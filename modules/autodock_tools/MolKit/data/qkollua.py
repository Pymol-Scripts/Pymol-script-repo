###################################################################
#  Kollman united atom (polar-Hydrogen) charges.
#
#  Garrett M. Morris, July 6, 1992. Original version in "awk".
#  Revised:        September 18, 1992,
#  Revised:        March 25, 1993.
#  Pythonised: Sometime in the 21st Century.
#
#
#  Now recognizes the alternate location indicator and 
#  chain-identifier correctly.
#

#
#  VERSION 3.0
#

#  Create a dictionary of partial atomic charges based on the 
#  Kollman United Atom Charge Set (KOLLUA), whose keys
#  are a tuple of the residue name and atom name:
#######################################################################

# $Header: /opt/cvs/python/packages/share1.5/MolKit/data/qkollua.py,v 1.1 2003/11/11 23:57:19 rhuey Exp $
#
# $Id: qkollua.py,v 1.1 2003/11/11 23:57:19 rhuey Exp $
#

q = {}

#
# ALA 6 0.000 Alanine, polypeptide residue
#
q[ "ALA"]={}
q[ "ALA"][ "N" ] = -0.5200
q[ "ALA"][ "CA" ] = 0.2150
q[ "ALA"][ "HN" ] = 0.2480
q[ "ALA"][ "C" ] = 0.5260
q[ "ALA"][ "O" ] = -0.5000
q[ "ALA"][ "CB" ] = 0.0310
#
# ARG 17 -12.480 Arginine(positive side chain), polypeptide residue
#
q[ "ARG"]={}
q[ "ARG"][ "N" ] = -0.5200
q[ "ARG"][ "CA" ] = 0.2370
q[ "ARG"][ "HN" ] = 0.2480
q[ "ARG"][ "C" ] = 0.5260
q[ "ARG"][ "O" ] = -0.5000
q[ "ARG"][ "CB" ] = 0.0490
q[ "ARG"][ "CG" ] = 0.0580
q[ "ARG"][ "CD" ] = 0.1110
q[ "ARG"][ "NE" ] = -0.4930
#adjusted charge on CZ by +.001
q[ "ARG"][ "CZ" ] = 0.8140
q[ "ARG"][ "NH1" ] = -0.6340
q[ "ARG"][ "NH2" ] = -0.6340
q[ "ARG"][ "HE" ] = 0.2940
q[ "ARG"][ "HH11" ] = 0.3610
q[ "ARG"][ "HH12" ] = 0.3610
q[ "ARG"][ "HH21" ] = 0.3610
q[ "ARG"][ "HH22" ] = 0.3610
#
# ASN 11 0.000 Asparagine, polypeptide residue
#
q[ "ASN"]={}
q[ "ASN"][ "N" ] = -0.5200
q[ "ASN"][ "CA" ] = 0.2170
q[ "ASN"][ "HN" ] = 0.2480
q[ "ASN"][ "C" ] = 0.5260
q[ "ASN"][ "O" ] = -0.5000
q[ "ASN"][ "CB" ] = 0.0030
q[ "ASN"][ "CG" ] = 0.6750
q[ "ASN"][ "OD1" ] = -0.4700
q[ "ASN"][ "ND2" ] = -0.8670
q[ "ASN"][ "HD21" ] = 0.3440
q[ "ASN"][ "HD22" ] = 0.3440
#
# ASPN 11 0.000 Aspartate(negative side chain), positive N-terminus
#
q[ "ASPN"]={}
q[ "ASPN"][ "N" ] = -0.2630
q[ "ASPN"][ "CA" ] = 0.3010
q[ "ASPN"][ "C" ] = 0.5260
q[ "ASPN"][ "O" ] = -0.5000
q[ "ASPN"][ "CB" ] = -0.2080
q[ "ASPN"][ "CG" ] = 0.6200
q[ "ASPN"][ "OD1" ] = -0.7060
q[ "ASPN"][ "OD2" ] = -0.7060
#
# ASP 9 3.650 Aspartate(negative side chain), polypeptide residue
#
q[ "ASP"]={}
q[ "ASP"][ "N" ] = -0.5200
q[ "ASP"][ "CA" ] = 0.2460
q[ "ASP"][ "HN" ] = 0.2480
q[ "ASP"][ "C" ] = 0.5260
q[ "ASP"][ "O" ] = -0.5000
q[ "ASP"][ "CB" ] = -0.2080
q[ "ASP"][ "CG" ] = 0.6200
q[ "ASP"][ "OD1" ] = -0.7060
q[ "ASP"][ "OD2" ] = -0.7060
#
# CYS 7 0.000 Cystine(sulfide bridge), polypeptide residue
#
q[ "CYS"]={}
q[ "CYS"][ "N" ] = -0.5200
q[ "CYS"][ "CA" ] = 0.0880
q[ "CYS"][ "HN" ] = 0.2480
q[ "CYS"][ "C" ] = 0.5260
q[ "CYS"][ "O" ] = -0.5000
q[ "CYS"][ "CB" ] = 0.1430
q[ "CYS"][ "SG" ] = 0.0150
#
# CYSH 8 0.000 Cysteine, polypeptide residue
#
q[ "CYSH"]={}
q[ "CYSH"][ "N" ] = -0.5200
q[ "CYSH"][ "CA" ] = 0.1460
q[ "CYSH"][ "HN" ] = 0.2480
q[ "CYSH"][ "C" ] = 0.5260
q[ "CYSH"][ "O" ] = -0.5000
q[ "CYSH"][ "CB" ] = 0.1000
q[ "CYSH"][ "SG" ] = -0.1350
q[ "CYSH"][ "HG" ] = 0.1350
#
# CYHC 9 0.000 Cysteine, negative C-terminus
#
q[ "CYHC"]={}
q[ "CYHC"][ "N" ] = -0.5200
q[ "CYHC"][ "CA" ] = 0.1400
q[ "CYHC"][ "HN" ] = 0.2480
q[ "CYHC"][ "C" ] = 0.4440
q[ "CYHC"][ "O" ] = -0.7060
q[ "CYHC"][ "CB" ] = 0.1000
q[ "CYHC"][ "SG" ] = -0.1350
q[ "CYHC"][ "HG" ] = 0.1350
#
# GLU 10 4.250 Glutamate(negative side chain), polypeptide residue
#
q[ "GLU"]={}
q[ "GLU"][ "N" ] = -0.5200
q[ "GLU"][ "CA" ] = 0.2460
q[ "GLU"][ "HN" ] = 0.2480
q[ "GLU"][ "C" ] = 0.5260
q[ "GLU"][ "O" ] = -0.5000
q[ "GLU"][ "CB" ] = -0.0000
q[ "GLU"][ "CG" ] = -0.2080
q[ "GLU"][ "CD" ] = 0.6200
q[ "GLU"][ "OE1" ] = -0.7060
q[ "GLU"][ "OE2" ] = -0.7060
#
# GLN 12 0.000 Glutamine, polypeptide residue
#
q[ "GLN"]={}
q[ "GLN"][ "N" ] = -0.5200
q[ "GLN"][ "CA" ] = 0.2100
q[ "GLN"][ "HN" ] = 0.2480
q[ "GLN"][ "C" ] = 0.5260
q[ "GLN"][ "O" ] = -0.5000
q[ "GLN"][ "CB" ] = 0.0530
q[ "GLN"][ "CG" ] = -0.0430
q[ "GLN"][ "CD" ] = 0.6750
q[ "GLN"][ "OE1" ] = -0.4700
q[ "GLN"][ "NE2" ] = -0.8670
q[ "GLN"][ "HE21" ] = 0.3440
q[ "GLN"][ "HE22" ] = 0.3440
#
# GLY 5 0.000 Glycine, polypeptide residue
#
q[ "GLY"]={}
q[ "GLY"][ "N" ] = -0.5200
q[ "GLY"][ "CA" ] = 0.2460
q[ "GLY"][ "HN" ] = 0.2480
q[ "GLY"][ "C" ] = 0.5260
q[ "GLY"][ "O" ] = -0.5000
#
# HIS 12 6.000 Histidine(H on NE),polypeptide residue
#
q[ "HIS"]={}
q[ "HIS"][ "N" ] = -0.5200
q[ "HIS"][ "CA" ] = 0.2190
q[ "HIS"][ "HN" ] = 0.2480
q[ "HIS"][ "C" ] = 0.5260
q[ "HIS"][ "O" ] = -0.5000
q[ "HIS"][ "CB" ] = 0.0600
q[ "HIS"][ "CG" ] = 0.1120
q[ "HIS"][ "ND1" ] = -0.5270
q[ "HIS"][ "CE1" ] = 0.3840
q[ "HIS"][ "NE2" ] = -0.4440
q[ "HIS"][ "CD2" ] = 0.1220
q[ "HIS"][ "HE2" ] = 0.3200
# other name variants: amber, lpdb
q[ "HIE"] = q[ "HSE"] = q[ "HIS"] 
#
# HIS+ 13 -6.000 Histidine(positive side chain), polypeptide residue
#
q[ "HIS+"]={}
q[ "HIS+"][ "N" ] = -0.5200
q[ "HIS+"][ "CA" ] = 0.1950
q[ "HIS+"][ "HN" ] = 0.2480
q[ "HIS+"][ "C" ] = 0.5260
q[ "HIS+"][ "O" ] = -0.5000
q[ "HIS+"][ "CB" ] = 0.2110
q[ "HIS+"][ "CG" ] = 0.1030
q[ "HIS+"][ "ND1" ] = -0.6130
q[ "HIS+"][ "CE1" ] = 0.7190
q[ "HIS+"][ "NE2" ] = -0.6860
q[ "HIS+"][ "CD2" ] = 0.3530
q[ "HIS+"][ "HD1" ] = 0.4780
q[ "HIS+"][ "HE2" ] = 0.4860
# other name variants: amber, lpdb
q[ "HIP"] = q[ "HSP"] = q[ "HIS+"] 
#
# HISD 12 6.000 Histidine(H on ND),polypeptide residue
#
q[ "HISD"]={}
q[ "HISD"][ "N" ] = -0.5200
q[ "HISD"][ "CA" ] = 0.2190
q[ "HISD"][ "HN" ] = 0.2480
q[ "HISD"][ "C" ] = 0.5260
q[ "HISD"][ "O" ] = -0.5000
q[ "HISD"][ "CB" ] = 0.0600
q[ "HISD"][ "CG" ] = 0.0890
q[ "HISD"][ "ND1" ] = -0.4440
q[ "HISD"][ "CE1" ] = 0.3840
q[ "HISD"][ "NE2" ] = -0.5270
q[ "HISD"][ "CD2" ] = 0.1450
q[ "HISD"][ "HD1" ] = 0.3200
# other name variants: amber, lpdb
q[ "HID"] = q[ "HSD"] = q[ "HISD"] 
#
# ILE 9 0.000 Isoleucine, polypeptide residue
#
q[ "ILE"]={}
q[ "ILE"][ "N" ] = -0.5200
q[ "ILE"][ "CA" ] = 0.1990
q[ "ILE"][ "HN" ] = 0.2480
q[ "ILE"][ "C" ] = 0.5260
q[ "ILE"][ "O" ] = -0.5000
q[ "ILE"][ "CB" ] = 0.0300
q[ "ILE"][ "CG2" ] = 0.0010
q[ "ILE"][ "CG1" ] = 0.0170
q[ "ILE"][ "CD1" ] = -0.0010
#
# ILEC 10 0.000 Isoleucine, negative C-terminus
#
q[ "ILEC"]={}
q[ "ILEC"][ "N" ] = -0.5200
q[ "ILEC"][ "CA" ] = 0.1930
q[ "ILEC"][ "HN" ] = 0.2480
q[ "ILEC"][ "C" ] = 0.4440
q[ "ILEC"][ "O" ] = -0.7060
q[ "ILEC"][ "CB" ] = 0.0300
q[ "ILEC"][ "CG2" ] = 0.0010
q[ "ILEC"][ "CG1" ] = 0.0170
q[ "ILEC"][ "CD1" ] = -0.0010
#
# LEU 9 0.000 Leucine, polypeptide residue
#
q[ "LEU"]={}
q[ "LEU"][ "N" ] = -0.5200
q[ "LEU"][ "CA" ] = 0.2040
q[ "LEU"][ "HN" ] = 0.2480
q[ "LEU"][ "C" ] = 0.5260
q[ "LEU"][ "O" ] = -0.5000
q[ "LEU"][ "CB" ] = 0.0160
q[ "LEU"][ "CG" ] = 0.0540
q[ "LEU"][ "CD1" ] = -0.0140
q[ "LEU"][ "CD2" ] = -0.0140
#
# LEUC 10 0.000 Leucine, negative C-terminus
#
q[ "LEUC"]={}
q[ "LEUC"][ "N" ] = -0.5200
q[ "LEUC"][ "CA" ] = 0.1980
q[ "LEUC"][ "HN" ] = 0.2480
q[ "LEUC"][ "C" ] = 0.4440
q[ "LEUC"][ "O" ] = -0.7060
q[ "LEUC"][ "CB" ] = 0.0160
q[ "LEUC"][ "CG" ] = 0.0540
q[ "LEUC"][ "CD1" ] = -0.0140
q[ "LEUC"][ "CD2" ] = -0.0140
#
# LYS 13 -10.530 Lysine(positive side chain), polypeptide residue
#
q[ "LYS"]={}
q[ "LYS"][ "N" ] = -0.5200
q[ "LYS"][ "CA" ] = 0.2270
q[ "LYS"][ "HN" ] = 0.2480
q[ "LYS"][ "C" ] = 0.5260
q[ "LYS"][ "O" ] = -0.5000
q[ "LYS"][ "CB" ] = 0.0390
q[ "LYS"][ "CG" ] = 0.0530
q[ "LYS"][ "CD" ] = 0.0480
q[ "LYS"][ "CE" ] = 0.2180
q[ "LYS"][ "NZ" ] = -0.2720
q[ "LYS"][ "HZ1" ] = 0.3110
q[ "LYS"][ "HZ2" ] = 0.3110
q[ "LYS"][ "HZ3" ] = 0.3110
#
# LYSN 15 0.000 Lysine(positive side chain), positive N-terminus
#
q[ "LYSN"]={}
q[ "LYSN"][ "N" ] = -0.2630
q[ "LYSN"][ "CA" ] = 0.2820
q[ "LYSN"][ "C" ] = 0.5260
q[ "LYSN"][ "O" ] = -0.5000
q[ "LYSN"][ "CB" ] = 0.0390
q[ "LYSN"][ "CG" ] = 0.0530
q[ "LYSN"][ "CD" ] = 0.0480
q[ "LYSN"][ "CE" ] = 0.2180
q[ "LYSN"][ "NZ" ] = -0.2720
q[ "LYSN"][ "HZ1" ] = 0.3110
q[ "LYSN"][ "HZ2" ] = 0.3110
q[ "LYSN"][ "HZ3" ] = 0.3110
#
# MET 9 0.000 Methionine, polypeptide residue
#
q[ "MET"]={}
q[ "MET"][ "N" ] = -0.5200
q[ "MET"][ "CA" ] = 0.1370
q[ "MET"][ "HN" ] = 0.2480
q[ "MET"][ "C" ] = 0.5260
q[ "MET"][ "O" ] = -0.5000
q[ "MET"][ "CB" ] = 0.0370
q[ "MET"][ "CG" ] = 0.0900
q[ "MET"][ "SD" ] = -0.0250
q[ "MET"][ "CE" ] = 0.0070
#
# PCA 9 0.000 cyclized gln, polypeptide residue
#
q[ "PCA"]={}
q[ "PCA"][ "N" ] = -0.5200
q[ "PCA"][ "CA" ] = 0.2100
q[ "PCA"][ "CB" ] = 0.0530
q[ "PCA"][ "CG" ] = -0.0430
q[ "PCA"][ "CD" ] = 0.5260
q[ "PCA"][ "OE" ] = -0.5000
q[ "PCA"][ "HN" ] = 0.2480
q[ "PCA"][ "C" ] = 0.5260
q[ "PCA"][ "O" ] = -0.5000
#
# PHE 12 0.000 Phenylalanine, polypeptide residue
#
q[ "PHE"]={}
q[ "PHE"][ "N" ] = -0.5200
q[ "PHE"][ "CA" ] = 0.2140
q[ "PHE"][ "HN" ] = 0.2480
q[ "PHE"][ "C" ] = 0.5260
q[ "PHE"][ "O" ] = -0.5000
q[ "PHE"][ "CB" ] = 0.0380
q[ "PHE"][ "CG" ] = 0.0110
q[ "PHE"][ "CD1" ] = -0.0110
q[ "PHE"][ "CE1" ] = 0.0040
q[ "PHE"][ "CZ" ] = -0.0030
q[ "PHE"][ "CE2" ] = 0.0040
q[ "PHE"][ "CD2" ] = -0.0110
#
# PRO 7 0.000 Proline, polypeptide residue
#
q[ "PRO"]={}
q[ "PRO"][ "N" ] = -0.2570
q[ "PRO"][ "CA" ] = 0.1120
q[ "PRO"][ "CD" ] = 0.0840
q[ "PRO"][ "C" ] = 0.5260
q[ "PRO"][ "O" ] = -0.5000
q[ "PRO"][ "CB" ] = -0.0010
q[ "PRO"][ "CG" ] = 0.0360
#
# SER 8 0.000 Serine, polypeptide residue
#
q[ "SER"]={}
q[ "SER"][ "N" ] = -0.5200
q[ "SER"][ "CA" ] = 0.2920
q[ "SER"][ "HN" ] = 0.2480
q[ "SER"][ "C" ] = 0.5260
q[ "SER"][ "O" ] = -0.5000
q[ "SER"][ "CB" ] = 0.1940
q[ "SER"][ "OG" ] = -0.5500
q[ "SER"][ "HG" ] = 0.3100
#
# THR 9 0.000 Threonine, polypeptide residue
#
q[ "THR"]={}
q[ "THR"][ "N" ] = -0.5200
q[ "THR"][ "CA" ] = 0.2680
q[ "THR"][ "HN" ] = 0.2480
q[ "THR"][ "C" ] = 0.5260
q[ "THR"][ "O" ] = -0.5000
q[ "THR"][ "CB" ] = 0.2110
q[ "THR"][ "OG1" ] = -0.5500
q[ "THR"][ "HG1" ] = 0.3100
q[ "THR"][ "CG2" ] = 0.0070
#
# TRP 16 0.000 Tryptophan, polypeptide residue
#
q[ "TRP"]={}
q[ "TRP"][ "N" ] = -0.5200
q[ "TRP"][ "CA" ] = 0.2480
q[ "TRP"][ "HN" ] = 0.2480
q[ "TRP"][ "C" ] = 0.5260
q[ "TRP"][ "O" ] = -0.5000
q[ "TRP"][ "CB" ] = 0.0200
q[ "TRP"][ "CG" ] = 0.0460
q[ "TRP"][ "CD1" ] = 0.1170
q[ "TRP"][ "NE1" ] = -0.3300
q[ "TRP"][ "CE2" ] = 0.0000
q[ "TRP"][ "CD2" ] = -0.2750
q[ "TRP"][ "HE1" ] = 0.2940
q[ "TRP"][ "CE3" ] = 0.1450
q[ "TRP"][ "CZ3" ] = -0.0820
q[ "TRP"][ "CH2" ] = 0.0340
q[ "TRP"][ "CZ2" ] = 0.0290
#
# TYR 14 -10.070 Tyrosine, polypeptide residue
#
q[ "TYR"]={}
q[ "TYR"][ "N" ] = -0.5200
q[ "TYR"][ "CA" ] = 0.2450
q[ "TYR"][ "HN" ] = 0.2480
q[ "TYR"][ "C" ] = 0.5260
q[ "TYR"][ "O" ] = -0.5000
q[ "TYR"][ "CB" ] = 0.0220
q[ "TYR"][ "CG" ] = -0.0010
q[ "TYR"][ "CD1" ] = -0.0350
q[ "TYR"][ "CE1" ] = 0.1000
q[ "TYR"][ "CZ" ] = -0.1210
q[ "TYR"][ "OH" ] = -0.3680
q[ "TYR"][ "HH" ] = 0.3390
q[ "TYR"][ "CE2" ] = 0.1000
q[ "TYR"][ "CD2" ] = -0.0350
#
# VAL 8 0.000 Valine, polypeptide residue
#
q[ "VAL"]={}
q[ "VAL"][ "N" ] = -0.5200
q[ "VAL"][ "CA" ] = 0.2010
q[ "VAL"][ "HN" ] = 0.2480
q[ "VAL"][ "C" ] = 0.5260
q[ "VAL"][ "O" ] = -0.5000
q[ "VAL"][ "CB" ] = 0.0330
q[ "VAL"][ "CG1" ] = 0.0060
q[ "VAL"][ "CG2" ] = 0.0060
#
# ACE 3 0.000 Acetyl, on N-terminus
#
q[ "ACE"]={}
q[ "ACE"][ "CA" ] = -0.0260
q[ "ACE"][ "C" ] = 0.5260
q[ "ACE"][ "O" ] = -0.5000
#
# N-M 3 0.000 N-Methyl, on C-terminus
#
q[ "N-M"]={}
q[ "N-M"][ "N" ] = -0.5000
q[ "N-M"][ "CA" ] = 0.2200
q[ "N-M"][ "HN" ] = 0.2800
#
# NH2 3 0.000 NH2, on C-terminus
#
q[ "NH2"]={}
q[ "NH2"][ "N" ] = -0.5000
q[ "NH2"][ "HN1" ] = 0.2500
q[ "NH2"][ "HN2" ] = 0.2500
#
# WTR 3 0.000 Water
#
q[ "WTR"]={}
q[ "WTR"][ "O" ] = -0.8200
q[ "WTR"][ "H1" ] = 0.4100
q[ "WTR"][ "H2" ] = 0.4100
q[ "WTR"][ "HO" ] = 0.4100
q[ "WTR"][ "HO1" ] = 0.4100
q[ "WTR"][ "HO2" ] = 0.4100
q[ "WTR"][ "O1" ] = -0.8200
q[ "WTR"][ "H11" ] = 0.4100
q[ "WTR"][ "H12" ] = 0.4100
#
# WTX 3 0.000 Water
#
q[ "WTX"]={}
q[ "WTX"][ "O" ] = -0.8200
q[ "WTX"][ "H1" ] = 0.4100
q[ "WTX"][ "H2" ] = 0.4100
q[ "WTX"][ "HO" ] = 0.4100
q[ "WTX"][ "HO1" ] = 0.4100
q[ "WTX"][ "HO2" ] = 0.4100
q[ "WTX"][ "O1" ] = -0.8200
q[ "WTX"][ "H11" ] = 0.4100
q[ "WTX"][ "H12" ] = 0.4100
#
# HOH 3 0.000 Water (RCSB Convention)
#
q[ "HOH"]={}
q[ "HOH"][ "O" ] = -0.8200
q[ "HOH"][ "H1" ] = 0.4100
q[ "HOH"][ "H2" ] = 0.4100
q[ "HOH"][ "HO" ] = 0.4100
q[ "HOH"][ "HO1" ] = 0.4100
q[ "HOH"][ "HO2" ] = 0.4100
q[ "HOH"][ "O1" ] = -0.8200
q[ "HOH"][ "H11" ] = 0.4100
q[ "HOH"][ "H12" ] = 0.4100
#
# WAT 3 0.000 Water (RCSB Convention)
#
q[ "WAT"]={}
q[ "WAT"][ "O" ] = -0.8200
q[ "WAT"][ "H1" ] = 0.4100
q[ "WAT"][ "H2" ] = 0.4100
q[ "WAT"][ "HO" ] = 0.4100
q[ "WAT"][ "HO1" ] = 0.4100
q[ "WAT"][ "HO2" ] = 0.4100
q[ "WAT"][ "O1" ] = -0.8200
q[ "WAT"][ "H11" ] = 0.4100
q[ "WAT"][ "H12" ] = 0.4100
#
# elib
#
#
#the following charges are from uni.in, part of the amber data set
#NB: it names the oxygen in the sugar ring 'O1*' whereas sol_par names it 'O4*'
#  D-ADENOSINE
q["DADE"] = q["  A"] = {}
q["DADE"]["O5*"]= -0.53500
q["DADE"]["O5'"]= -0.53500
q["DADE"]["C5*"]=  0.15300
q["DADE"]["C5'"]=  0.15300
q["DADE"]["C4*"]=  0.18500
q["DADE"]["C4'"]=  0.18500
q["DADE"]["O1*"]= -0.38600
q["DADE"]["O1'"]= -0.38600
q["DADE"]["O4*"]= -0.38600
q["DADE"]["O4'"]= -0.38600
q["DADE"]["C1*"]=  0.50000
q["DADE"]["C1'"]=  0.50000
q["DADE"]["N9"]= -0.45700
q["DADE"]["C8"]=  0.48800
q["DADE"]["N7"]= -0.59900
q["DADE"]["C5"]= -0.15100
q["DADE"]["C6"]=  0.81300
q["DADE"]["N6"]= -0.79300
q["DADE"]["H61"]=  0.33500
q["DADE"]["H62"]=  0.33900
q["DADE"]["N1"]= -0.76000
q["DADE"]["C2"]=  0.57100
q["DADE"]["N3"]= -0.71700
q["DADE"]["C4"]=  0.69500
q["DADE"]["C3*"]=  0.17200
q["DADE"]["C3'"]=  0.17200
q["DADE"]["C2*"]= -0.04700
q["DADE"]["C2'"]= -0.04700
q["DADE"]["O3*"]= -0.53500
q["DADE"]["O3'"]= -0.53500
q["DADE"]["P"] =   1.42900
q["DADE"]["O1P"] =  -0.8500
q["DADE"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["DADE"]["H1P"] =  0.00
q["DADE"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["DADE"]["H5*"] =  0.00
q["DADE"]["H5'"] =  0.00
q["DADE"]["H5*1"] =  0.00
q["DADE"]["H5'1"] =  0.00
q["DADE"]["H5*2"] =  0.00
q["DADE"]["H5'2"] =  0.00
q["DADE"]["H4*"] =  0.00
q["DADE"]["H4'"] =  0.00
q["DADE"]["H3*"] =  0.00
q["DADE"]["H3'"] =  0.00
q["DADE"]["H2*1"] =  0.00
q["DADE"]["H2'1"] =  0.00
q["DADE"]["H2*2"] =  0.00
q["DADE"]["H2'2"] =  0.00
q["DADE"]["H1*"] =  0.00
q["DADE"]["H1'"] =  0.00
q["DADE"]["H8"] =  0.00
q["DADE"]["H7"] =  0.00
q["DADE"]["H3"] =  0.00
q["DADE"]["H2"] =  0.00
q["DADE"]["H1"] =  0.00

#  R-ADENOSINE
q["RADE"] = q[" RA"] = {}
q["RADE"]["O5*"] = -0.53500
q["RADE"]["O5'"] = -0.53500
q["RADE"]["C5*"] =  0.17500
q["RADE"]["C5'"] =  0.17500
q["RADE"]["C4*"] =  0.19700
q["RADE"]["C4'"] =  0.19700
q["RADE"]["O1*"] = -0.41300
q["RADE"]["O1'"] = -0.41300
q["RADE"]["O4*"] = -0.41300
q["RADE"]["O4'"] = -0.41300
q["RADE"]["C1*"] =  0.52200
q["RADE"]["C1'"] =  0.52200
q["RADE"]["N9"]  = -0.45700
q["RADE"]["C8"]  =  0.48800
q["RADE"]["N7"]  = -0.59900
q["RADE"]["C5"]  = -0.15100
q["RADE"]["C6"]  =  0.81300
q["RADE"]["N6"]  = -0.79300
q["RADE"]["H61"] =  0.33500
q["RADE"]["H62"] =  0.33900
q["RADE"]["N1"]  = -0.76000
q["RADE"]["C2"]  =  0.57100
q["RADE"]["N3"]  = -0.71700
q["RADE"]["C4"]  =  0.69500
q["RADE"]["C3*"] =  0.21000
q["RADE"]["C3'"] =  0.21000
q["RADE"]["C2*"] =  0.08200
q["RADE"]["C2'"] =  0.08200
q["RADE"]["O2*"] = -0.51200    #for rna only
q["RADE"]["O2'"] = -0.51200    #for rna only
q["RADE"]["HO2*"] =  0.31600    #for rna only
q["RADE"]["HO2'"] =  0.31600    #for rna only
q["RADE"]["O3*"] = -0.53500
q["RADE"]["O3'"] = -0.53500
q["RADE"]["P"] =   1.42900
q["RADE"]["O1P"] =  -0.8500
q["RADE"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["RADE"]["H1P"] =  0.00
q["RADE"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["RADE"]["H5*"] =  0.00
q["RADE"]["H5'"] =  0.00
q["RADE"]["H5*1"] =  0.00
q["RADE"]["H5'1"] =  0.00
q["RADE"]["H5*2"] =  0.00
q["RADE"]["H5'2"] =  0.00
q["RADE"]["H4*"] =  0.00
q["RADE"]["H4'"] =  0.00
q["RADE"]["H3*"] =  0.00
q["RADE"]["H3'"] =  0.00
#if O2* atom, C2* could still have 1 Hydrogen
q["RADE"]["H2*"] =  0.00
q["RADE"]["H2'"] =  0.00
q["RADE"]["H1*"] =  0.00
q["RADE"]["H1'"] =  0.00
q["RADE"]["H8"] =  0.00
q["RADE"]["H7"] =  0.00
q["RADE"]["H3"] =  0.00
q["RADE"]["H2"] =  0.00
q["RADE"]["H1"] =  0.00

#  D-THYMINE
q["DTHY"] = q["  T"] = {}
q["DTHY"]["O5*"] =-0.53500
q["DTHY"]["O5'"] =-0.53500
q["DTHY"]["C5*"] = 0.15300
q["DTHY"]["C5'"] = 0.15300
q["DTHY"]["C4*"] = 0.18500
q["DTHY"]["C4'"] = 0.18500
q["DTHY"]["O1*"] =-0.38600
q["DTHY"]["O1'"] =-0.38600
q["DTHY"]["O4*"] =-0.38600
q["DTHY"]["O4'"] =-0.38600
q["DTHY"]["C1*"] = 0.50000
q["DTHY"]["C1'"] = 0.50000
q["DTHY"]["N1"] = -0.73900
q["DTHY"]["C6"] =  0.55100
q["DTHY"]["C5"] = -0.59500
q["DTHY"]["C7"] =  0.09700
q["DTHY"]["C5M"] =  0.09700
q["DTHY"]["C4"] =  0.98000
q["DTHY"]["O4"] = -0.47200
q["DTHY"]["N3"] = -1.01200
q["DTHY"]["H3"] =  0.37000
q["DTHY"]["C2"] =  1.11300
q["DTHY"]["O2"] = -0.52900
q["DTHY"]["C3*"] = 0.17200
q["DTHY"]["C3'"] = 0.17200
q["DTHY"]["C2*"] =-0.04700
q["DTHY"]["C2'"] =-0.04700
q["DTHY"]["O3*"] =-0.53500
q["DTHY"]["O3'"] =-0.53500
q["DTHY"]["P"] =   1.42900
q["DTHY"]["O1P"] =  -0.8500
q["DTHY"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["DTHY"]["H1P"] =  0.00
q["DTHY"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["DTHY"]["H5*1"] =  0.00
q["DTHY"]["H5'1"] =  0.00
q["DTHY"]["H5*2"] =  0.00
q["DTHY"]["H5'2"] =  0.00
q["DTHY"]["H4*"] =  0.00
q["DTHY"]["H4'"] =  0.00
q["DTHY"]["H3*"] =  0.00
q["DTHY"]["H3'"] =  0.00
q["DTHY"]["H2*1"] =  0.00
q["DTHY"]["H2'1"] =  0.00
q["DTHY"]["H2*2"] =  0.00
q["DTHY"]["H2'2"] =  0.00
q["DTHY"]["H1*"] =  0.00
q["DTHY"]["H1'"] =  0.00
q["DTHY"]["H6"] =  0.00
q["DTHY"]["H5M1"] =  0.00
q["DTHY"]["H5M2"] =  0.00
q["DTHY"]["H5M3"] =  0.00

#  R-URACIL
q["RURA"] = q["  U"] = {}
q["RURA"]["O5*"] =  -0.53500
q["RURA"]["O5'"] =  -0.53500
q["RURA"]["C5*"] =   0.17500
q["RURA"]["C5'"] =   0.17500
q["RURA"]["C4*"] =   0.19700
q["RURA"]["C4'"] =   0.19700
q["RURA"]["O1*"] =  -0.41300
q["RURA"]["O1'"] =  -0.41300
q["RURA"]["O4*"] =  -0.41300
q["RURA"]["O4'"] =  -0.41300
q["RURA"]["C1*"] =   0.52200
q["RURA"]["C1'"] =   0.52200
q["RURA"]["N1"] =   -0.56700
q["RURA"]["C6"] =    0.36600
q["RURA"]["C5"] =   -0.20400
q["RURA"]["C4"] =    0.57200
q["RURA"]["O4"] =   -0.39400
q["RURA"]["N3"] =   -0.75800
q["RURA"]["H3"] =    0.34700
q["RURA"]["C2"] =    0.89600
q["RURA"]["O2"] =   -0.49400
q["RURA"]["C3*"] =   0.21000
q["RURA"]["C3'"] =   0.21000
q["RURA"]["C2*"] =   0.08200
q["RURA"]["C2'"] =   0.08200
q["RURA"]["O2*"] =  -0.51200
q["RURA"]["O2'"] =  -0.51200
q["RURA"]["HO2*"] =  0.31600
q["RURA"]["HO2'"] =  0.31600
q["RURA"]["O3*"] =  -0.53500
q["RURA"]["O3'"] =  -0.53500
q["RURA"]["P"] =   1.42900
q["RURA"]["O1P"] =  -0.8500
q["RURA"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["RURA"]["H1P"] =  0.00
q["RURA"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
#the sugar's nphs
q["RURA"]["H5*1"] =  0.00
q["RURA"]["H5'1"] =  0.00
q["RURA"]["H5*2"] =  0.00
q["RURA"]["H5'2"] =  0.00
q["RURA"]["H4*"] =  0.00
q["RURA"]["H4'"] =  0.00
q["RURA"]["H3*"] =  0.00
q["RURA"]["H3'"] =  0.00
#if O2* atom, C2* could still have 1 Hydrogen
q["RURA"]["H2*"] =  0.00
q["RURA"]["H2'"] =  0.00
q["RURA"]["H1*"] =  0.00
q["RURA"]["H1'"] =  0.00
#the base's nphs
q["RURA"]["H51"] =  0.00
q["RURA"]["H52"] =  0.00
q["RURA"]["H61"] =  0.00
q["RURA"]["H62"] =  0.00

#  D-GUANOSINE
q["DGUA"] = q["  G"] = {}
q["DGUA"]["O5*"] =  -0.53500
q["DGUA"]["O5'"] =  -0.53500
q["DGUA"]["C5*"] =   0.15300
q["DGUA"]["C5'"] =   0.15300
q["DGUA"]["C4*"] =   0.18500
q["DGUA"]["C4'"] =   0.18500
q["DGUA"]["O1*"] =  -0.38600
q["DGUA"]["O1'"] =  -0.38600
q["DGUA"]["O4*"] =  -0.38600
q["DGUA"]["O4'"] =  -0.38600
q["DGUA"]["C1*"] =   0.50000
q["DGUA"]["C1'"] =   0.50000
q["DGUA"]["N9"] =   -0.37900
q["DGUA"]["C8"] =    0.42800
q["DGUA"]["N7"] =   -0.57500
q["DGUA"]["C5"] =   -0.08800
q["DGUA"]["C6"] =    0.71400
q["DGUA"]["O6"] =   -0.45900
q["DGUA"]["N1"] =   -0.74600
q["DGUA"]["H1"] =    0.34000
q["DGUA"]["C2"] =    0.84200
q["DGUA"]["N2"] =   -0.75800
q["DGUA"]["H21"] =   0.32400
q["DGUA"]["H22"] =   0.33300
q["DGUA"]["N3"] =  -0.70200
q["DGUA"]["C4"] =    0.49000
q["DGUA"]["C3*"] =   0.17200
q["DGUA"]["C3'"] =   0.17200
q["DGUA"]["C2*"] =  -0.04700
q["DGUA"]["C2'"] =  -0.04700
q["DGUA"]["O3*"] =  -0.53500
q["DGUA"]["O3'"] =  -0.53500
q["DGUA"]["P"] =   1.42900
q["DGUA"]["O1P"] =  -0.8500
q["DGUA"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["DGUA"]["H1P"] =  0.00
q["DGUA"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["DGUA"]["H5*1"] =  0.00
q["DGUA"]["H5'1"] =  0.00
q["DGUA"]["H5*2"] =  0.00
q["DGUA"]["H5'2"] =  0.00
q["DGUA"]["H4*"] =  0.00
q["DGUA"]["H4'"] =  0.00
q["DGUA"]["H3*"] =  0.00
q["DGUA"]["H3'"] =  0.00
q["DGUA"]["H2*1"] =  0.00
q["DGUA"]["H2'1"] =  0.00
q["DGUA"]["H2*2"] =  0.00
q["DGUA"]["H2'2"] =  0.00
q["DGUA"]["H1*"] =  0.00
q["DGUA"]["H1'"] =  0.00
q["DGUA"]["H8"] =  0.00

#  R-GUANOSINE
q["RGUA"] = q[" RG"] = {}
q["RGUA"]["O5*"] =  -0.53500
q["RGUA"]["O5'"] =  -0.53500
q["RGUA"]["C5*"] =   0.17500
q["RGUA"]["C5'"] =   0.17500
q["RGUA"]["C4*"] =   0.19700
q["RGUA"]["C4'"] =   0.19700
q["RGUA"]["O1*"] =  -0.41300
q["RGUA"]["O1'"] =  -0.41300
q["RGUA"]["O4*"] =  -0.41300
q["RGUA"]["O4'"] =  -0.41300
q["RGUA"]["C1*"] =   0.52200
q["RGUA"]["C1'"] =   0.52200
q["RGUA"]["N9"] =   -0.37900
q["RGUA"]["C8"] =    0.42800
q["RGUA"]["N7"] =   -0.57500
q["RGUA"]["C5"] =   -0.08800
q["RGUA"]["C6"] =    0.71400
q["RGUA"]["O6"] =   -0.45900
q["RGUA"]["N1"] =   -0.74600
q["RGUA"]["H1"] =    0.34000
q["RGUA"]["C2"] =    0.84200
q["RGUA"]["N2"] =   -0.75800
q["RGUA"]["H21"] =   0.32400
q["RGUA"]["H22"] =   0.33300
q["RGUA"]["N3"] =   -0.70200
q["RGUA"]["C4"] =    0.49000
q["RGUA"]["C3*"] =   0.21000
q["RGUA"]["C3'"] =   0.21000
q["RGUA"]["C2*"] =   0.08200
q["RGUA"]["C2'"] =   0.08200
q["RGUA"]["O2*"] =  -0.51200
q["RGUA"]["O2'"] =  -0.51200
q["RGUA"]["HO2*"] =  0.31600
q["RGUA"]["HO2'"] =  0.31600
q["RGUA"]["O3*"] =  -0.53500
q["RGUA"]["O3'"] =  -0.53500
q["RGUA"]["P"] =   1.42900
q["RGUA"]["O1P"] =  -0.8500
q["RGUA"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["RGUA"]["H1P"] =  0.00
q["RGUA"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["RGUA"]["H5*1"] =  0.00
q["RGUA"]["H5'1"] =  0.00
q["RGUA"]["H5*2"] =  0.00
q["RGUA"]["H5'2"] =  0.00
q["RGUA"]["H4*"] =  0.00
q["RGUA"]["H4'"] =  0.00
q["RGUA"]["H3*"] =  0.00
q["RGUA"]["H3'"] =  0.00
#if O2* atom, C2* could still have H2*
q["RGUA"]["H2*"] =  0.00
q["RGUA"]["H2'"] =  0.00
q["RGUA"]["H1*"] =  0.00
q["RGUA"]["H1'"] =  0.00
q["RGUA"]["H8"] =  0.00

#  D-CYTOSINE
q["DCYT"] = q["  C"] = {}
q["DCYT"]["O5*"] =  -0.53500
q["DCYT"]["O5'"] =  -0.53500
q["DCYT"]["C5*"] =   0.15300
q["DCYT"]["C5'"] =   0.15300
q["DCYT"]["C4*"] =   0.18500
q["DCYT"]["C4'"] =   0.18500
q["DCYT"]["O1*"] =  -0.38600
q["DCYT"]["O1'"] =  -0.38600
q["DCYT"]["O4*"] =  -0.38600
q["DCYT"]["O4'"] =  -0.38600
q["DCYT"]["C1*"] =   0.50000
q["DCYT"]["C1'"] =   0.50000
q["DCYT"]["N1"] =   -0.57200
q["DCYT"]["C6"] =    0.37700
q["DCYT"]["C5"] =   -0.23000
q["DCYT"]["C4"] =    0.63000
q["DCYT"]["N4"] =   -0.74300
q["DCYT"]["H41"] =   0.33500
q["DCYT"]["H42"] =   0.33800
q["DCYT"]["N3"] =   -0.79100
q["DCYT"]["C2"] =    0.93800
q["DCYT"]["O2"] =   -0.51800
q["DCYT"]["C3*"] =   0.17200
q["DCYT"]["C3'"] =   0.17200
q["DCYT"]["C2*"] =  -0.04700
q["DCYT"]["C2'"] =  -0.04700
q["DCYT"]["O3*"] =  -0.53500
q["DCYT"]["O3'"] =  -0.53500
q["DCYT"]["P"] =   1.42900
q["DCYT"]["O1P"] =  -0.8500
q["DCYT"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["DCYT"]["H1P"] =  0.00
q["DCYT"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["DCYT"]["H5*1"] =  0.00
q["DCYT"]["H5'1"] =  0.00
q["DCYT"]["H5*2"] =  0.00
q["DCYT"]["H5'2"] =  0.00
q["DCYT"]["H4*"] =  0.00
q["DCYT"]["H4'"] =  0.00
q["DCYT"]["H3*"] =  0.00
q["DCYT"]["H3'"] =  0.00
q["DCYT"]["H2*1"] =  0.00
q["DCYT"]["H2'1"] =  0.00
q["DCYT"]["H2*2"] =  0.00
q["DCYT"]["H2'2"] =  0.00
q["DCYT"]["H1*"] =  0.00
q["DCYT"]["H1'"] =  0.00
q["DCYT"]["H5"] =  0.00
q["DCYT"]["H6"] =  0.00

#   R-CYTOSINE
q["RCYT"] = q[" RC"] = {}
q["RCYT"]["O5*"] =  -0.53500
q["RCYT"]["O5'"] =  -0.53500
q["RCYT"]["C5*"] =   0.17500
q["RCYT"]["C5'"] =   0.17500
q["RCYT"]["C4*"] =   0.19700
q["RCYT"]["C4'"] =   0.19700
q["RCYT"]["O1*"] =  -0.41300
q["RCYT"]["O1'"] =  -0.41300
q["RCYT"]["O4*"] =  -0.41300
q["RCYT"]["O4'"] =  -0.41300
q["RCYT"]["C1*"] =   0.52200
q["RCYT"]["C1'"] =   0.52200
q["RCYT"]["N1"] =  -0.57200
q["RCYT"]["C6"] =    0.37700
q["RCYT"]["C5"] =   -0.23000
q["RCYT"]["C4"] =    0.63000
q["RCYT"]["N4"] =   -0.74300
q["RCYT"]["H41"] =   0.33500
q["RCYT"]["H42"] =   0.33800
q["RCYT"]["N3"] =   -0.79100
q["RCYT"]["C2"] =    0.93800
q["RCYT"]["O2"] =   -0.51800
q["RCYT"]["C3*"] =   0.21000
q["RCYT"]["C3'"] =   0.21000
q["RCYT"]["C2*"] =   0.08200
q["RCYT"]["C2'"] =   0.08200
q["RCYT"]["O2*"] =  -0.51200
q["RCYT"]["O2'"] =  -0.51200
q["RCYT"]["HO2*"] =  0.31600
q["RCYT"]["HO2'"] =  0.31600
q["RCYT"]["O3*"] =  -0.53500
q["RCYT"]["O3'"] =  -0.53500
q["RCYT"]["P"] =   1.42900
q["RCYT"]["O1P"] =  -0.8500
q["RCYT"]["O2P"] =  -0.8500
#no charges available for these hydrogens:
q["RCYT"]["H1P"] =  0.00
q["RCYT"]["H2P"] =  0.00
# possible other nph names generated by add_hydrogens
q["RCYT"]["H5*1"] =  0.00
q["RCYT"]["H5'1"] =  0.00
q["RCYT"]["H5*2"] =  0.00
q["RCYT"]["H5'2"] =  0.00
q["RCYT"]["H4*"] =  0.00
q["RCYT"]["H4'"] =  0.00
q["RCYT"]["H3*"] =  0.00
q["RCYT"]["H3'"] =  0.00
#if O2* atom, C2* could still have 1 Hydrogen
q["RCYT"]["H2*"] =  0.00
q["RCYT"]["H2'"] =  0.00
q["RCYT"]["H1*"] =  0.00
q["RCYT"]["H1'"] =  0.00
q["RCYT"]["H5"] =  0.00
q["RCYT"]["H6"] =  0.00

#  D - PHOSPATE MINUS
q["DPOM"] = {}
q["DPOM"]["P"] =  1.42900
q["DPOM"]["OA"] =-0.85000
q["DPOM"]["OB"] =-0.85000
# possible other nph names generated by add_hydrogens
q["DPOM"]["HOA"] =  0.00
q["DPOM"]["HOB"] =  0.00
q["DPOM"]["HOL"] =  0.00
q["DPOM"]["HOR"] =  0.00


#  R - PHOSPHATE MINUS
q["RPOM"] = {}
q["RPOM"]["P"] =  1.42900
q["RPOM"]["OA"] =-0.85000
q["RPOM"]["OB"] =-0.85000
# possible other nph names generated by add_hydrogens
q["RPOM"]["HOA"] =  0.00
q["RPOM"]["HOB"] =  0.00
q["RPOM"]["HOL"] =  0.00
q["RPOM"]["HOR"] =  0.00
