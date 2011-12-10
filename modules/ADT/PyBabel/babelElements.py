#############################################################################
#
# Author: Michel F. SANNER
#         Garrett M. MORRIS
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
# 2002/08/12 21:55:00 garrett:
# Elements sorted according to atomic number by GMM.
#
# 2002/08/12 21:55:00 garrett:
# guess_babel_elements added by GMM, to help make code more forgiving
# and less likely to croak.
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/babelElements.py,v 1.2 2002/09/11 18:15:15 garrett Exp $
#
# $Id: babelElements.py,v 1.2 2002/09/11 18:15:15 garrett Exp $
#


# File generated from element.lis
#
babel_elements = {
'Xx': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 0, 'rgb': (0.07, 0.5, 0.7)}, 
'H': {'bond_ord_rad': 0.352, 'cov_rad': 0.23, 'bs_rad': 0.24, 'vdw_rad': 1.2, 'max_bonds': 1, 'num': 1, 'rgb': (1.0, 1.0, 1.0)}, 
'He': {'bond_ord_rad': 0.7, 'cov_rad': 0.7, 'bs_rad': 0.244, 'vdw_rad': 1.22, 'max_bonds': 0, 'num': 2, 'rgb': (1.0, 0.37, 0.08)}, 
'Li': {'bond_ord_rad': 1.23, 'cov_rad': 0.68, 'bs_rad': 0.304, 'vdw_rad': 1.52, 'max_bonds': 1, 'num': 3, 'rgb': (0.0, 0.0, 0.1)}, 
'Be': {'bond_ord_rad': 0.9, 'cov_rad': 0.35, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 2, 'num': 4, 'rgb': (1.0, 0.37, 0.08)}, 
'B': {'bond_ord_rad': 1.0, 'cov_rad': 0.83, 'bs_rad': 0.416, 'vdw_rad': 2.08, 'max_bonds': 3, 'num': 5, 'rgb': (0.07, 0.5, 0.7)}, 
'C': {'bond_ord_rad': 0.762, 'cov_rad': 0.68, 'bs_rad': 0.37, 'vdw_rad': 1.85, 'max_bonds': 4, 'num': 6, 'rgb': (0.5, 0.5, 0.5)}, 
'N': {'bond_ord_rad': 0.676, 'cov_rad': 0.68, 'bs_rad': 0.308, 'vdw_rad': 1.54, 'max_bonds': 4, 'num': 7, 'rgb': (0.0, 0.0, 1.0)}, 
'O': {'bond_ord_rad': 0.64, 'cov_rad': 0.68, 'bs_rad': 0.28, 'vdw_rad': 1.4, 'max_bonds': 2, 'num': 8, 'rgb': (1.0, 0.0, 0.0)}, 
'F': {'bond_ord_rad': 0.63, 'cov_rad': 0.64, 'bs_rad': 0.27, 'vdw_rad': 1.35, 'max_bonds': 1, 'num': 9, 'rgb': (0.0, 0.0, 0.1)}, 
'Ne': {'bond_ord_rad': 0.7, 'cov_rad': 0.7, 'bs_rad': 0.32, 'vdw_rad': 1.6, 'max_bonds': 0, 'num': 10, 'rgb': (1.0, 0.37, 0.08)}, 
'Na': {'bond_ord_rad': 1.54, 'cov_rad': 0.97, 'bs_rad': 0.462, 'vdw_rad': 2.31, 'max_bonds': 1, 'num': 11, 'rgb': (0.0, 0.0, 0.1)}, 
'Mg': {'bond_ord_rad': 1.36, 'cov_rad': 1.1, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 2, 'num': 12, 'rgb': (0.0, 0.0, 0.1)}, 
'Al': {'bond_ord_rad': 1.18, 'cov_rad': 1.35, 'bs_rad': 0.41, 'vdw_rad': 2.05, 'max_bonds': 6, 'num': 13, 'rgb': (0.0, 0.0, 0.1)}, 
'Si': {'bond_ord_rad': 1.118, 'cov_rad': 1.2, 'bs_rad': 0.4, 'vdw_rad': 2.0, 'max_bonds': 6, 'num': 14, 'rgb': (0.7, 0.0, 0.1)}, 
'P': {'bond_ord_rad': 1.094, 'cov_rad': 0.75, 'bs_rad': 0.38, 'vdw_rad': 1.4, 'max_bonds': 5, 'num': 15, 'rgb': (0.5, 0.1, 0.8)}, 
'S': {'bond_ord_rad': 1.053, 'cov_rad': 1.02, 'bs_rad': 0.37, 'vdw_rad': 1.85, 'max_bonds': 6, 'num': 16, 'rgb': (1.0, 1.0, 0.0)}, 
'Cl': {'bond_ord_rad': 1.033, 'cov_rad': 0.99, 'bs_rad': 0.362, 'vdw_rad': 1.81, 'max_bonds': 1, 'num': 17, 'rgb': (0.0, 0.0, 0.1)}, 
'Ar': {'bond_ord_rad': 1.74, 'cov_rad': 0.7, 'bs_rad': 0.382, 'vdw_rad': 1.91, 'max_bonds': 0, 'num': 18, 'rgb': (1.0, 0.37, 0.08)}, 
'K': {'bond_ord_rad': 2.03, 'cov_rad': 1.33, 'bs_rad': 0.462, 'vdw_rad': 2.31, 'max_bonds': 1, 'num': 19, 'rgb': (0.0, 0.0, 0.1)}, 
'Ca': {'bond_ord_rad': 1.74, 'cov_rad': 0.99, 'bs_rad': 0.395, 'vdw_rad': 1.973, 'max_bonds': 2, 'num': 20, 'rgb': (0.0, 0.0, 0.1)}, 
'Sc': {'bond_ord_rad': 1.44, 'cov_rad': 1.44, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 21, 'rgb': (0.05, 0.05, 0.5)}, 
'Ti': {'bond_ord_rad': 1.32, 'cov_rad': 1.47, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 22, 'rgb': (0.05, 0.05, 0.5)}, 
'V': {'bond_ord_rad': 1.22, 'cov_rad': 1.33, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 23, 'rgb': (0.05, 0.05, 0.5)}, 
'Cr': {'bond_ord_rad': 1.18, 'cov_rad': 1.35, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 24, 'rgb': (0.05, 0.05, 0.5)}, 
'Mn': {'bond_ord_rad': 1.17, 'cov_rad': 1.35, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 8, 'num': 25, 'rgb': (0.05, 0.05, 0.5)}, 
'Fe': {'bond_ord_rad': 1.17, 'cov_rad': 1.34, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 26, 'rgb': (0.7, 0.0, 0.1)}, 
'Co': {'bond_ord_rad': 1.16, 'cov_rad': 1.33, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 27, 'rgb': (0.05, 0.05, 0.5)}, 
'Ni': {'bond_ord_rad': 1.15, 'cov_rad': 1.5, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 28, 'rgb': (0.05, 0.05, 0.5)}, 
'Cu': {'bond_ord_rad': 1.17, 'cov_rad': 1.52, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 29, 'rgb': (0.05, 0.05, 0.5)}, 
'Zn': {'bond_ord_rad': 1.25, 'cov_rad': 1.45, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 30, 'rgb': (0.05, 0.05, 0.5)}, 
'Ga': {'bond_ord_rad': 1.26, 'cov_rad': 1.22, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 3, 'num': 31, 'rgb': (0.05, 0.05, 0.5)}, 
'Ge': {'bond_ord_rad': 1.188, 'cov_rad': 1.17, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 4, 'num': 32, 'rgb': (0.05, 0.05, 0.5)}, 
'As': {'bond_ord_rad': 1.2, 'cov_rad': 1.21, 'bs_rad': 0.4, 'vdw_rad': 2.0, 'max_bonds': 3, 'num': 33, 'rgb': (0.4, 1.0, 0.1)}, 
'Se': {'bond_ord_rad': 1.187, 'cov_rad': 1.22, 'bs_rad': 0.4, 'vdw_rad': 2.0, 'max_bonds': 2, 'num': 34, 'rgb': (0.8, 1.0, 0.1)}, 
'Br': {'bond_ord_rad': 1.187, 'cov_rad': 1.21, 'bs_rad': 0.42, 'vdw_rad': 2.1, 'max_bonds': 1, 'num': 35, 'rgb': (0.0, 0.0, 0.1)}, 
'Kr': {'bond_ord_rad': 1.91, 'cov_rad': 1.91, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 0, 'num': 36, 'rgb': (1.0, 0.37, 0.08)}, 
'Rb': {'bond_ord_rad': 2.16, 'cov_rad': 1.47, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 1, 'num': 37, 'rgb': (0.05, 0.05, 0.5)}, 
'Sr': {'bond_ord_rad': 1.91, 'cov_rad': 1.12, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 2, 'num': 38, 'rgb': (0.05, 0.05, 0.5)}, 
'Y': {'bond_ord_rad': 1.62, 'cov_rad': 1.78, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 39, 'rgb': (0.05, 0.05, 0.5)}, 
'Zr': {'bond_ord_rad': 1.45, 'cov_rad': 1.56, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 40, 'rgb': (0.05, 0.05, 0.5)}, 
'Nb': {'bond_ord_rad': 1.34, 'cov_rad': 1.48, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 41, 'rgb': (0.05, 0.05, 0.5)}, 
'Mo': {'bond_ord_rad': 1.3, 'cov_rad': 1.47, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 42, 'rgb': (0.05, 0.05, 0.5)}, 
'Tc': {'bond_ord_rad': 1.27, 'cov_rad': 1.35, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 43, 'rgb': (0.05, 0.05, 0.5)}, 
'Ru': {'bond_ord_rad': 1.25, 'cov_rad': 1.4, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 44, 'rgb': (0.05, 0.05, 0.5)}, 
'Rh': {'bond_ord_rad': 1.25, 'cov_rad': 1.45, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 45, 'rgb': (0.05, 0.05, 0.5)}, 
'Pd': {'bond_ord_rad': 1.28, 'cov_rad': 1.5, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 46, 'rgb': (0.05, 0.05, 0.5)}, 
'Ag': {'bond_ord_rad': 1.34, 'cov_rad': 1.59, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 47, 'rgb': (0.05, 0.05, 0.5)}, 
'Cd': {'bond_ord_rad': 1.48, 'cov_rad': 1.69, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 48, 'rgb': (0.05, 0.05, 0.5)}, 
'In': {'bond_ord_rad': 1.44, 'cov_rad': 1.63, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 3, 'num': 49, 'rgb': (0.05, 0.05, 0.5)}, 
'Sn': {'bond_ord_rad': 1.385, 'cov_rad': 1.46, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 4, 'num': 50, 'rgb': (0.05, 0.05, 0.5)}, 
'Sb': {'bond_ord_rad': 1.4, 'cov_rad': 1.46, 'bs_rad': 0.44, 'vdw_rad': 2.2, 'max_bonds': 3, 'num': 51, 'rgb': (0.05, 0.05, 0.5)}, 
'Te': {'bond_ord_rad': 1.378, 'cov_rad': 1.47, 'bs_rad': 0.44, 'vdw_rad': 2.2, 'max_bonds': 2, 'num': 52, 'rgb': (0.05, 0.05, 0.5)}, 
'I': {'bond_ord_rad': 1.387, 'cov_rad': 1.4, 'bs_rad': 0.43, 'vdw_rad': 2.15, 'max_bonds': 1, 'num': 53, 'rgb': (0.05, 0.05, 0.5)}, 
'Xe': {'bond_ord_rad': 1.98, 'cov_rad': 1.98, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 0, 'num': 54, 'rgb': (0.05, 0.05, 0.5)}, 
'Cs': {'bond_ord_rad': 2.35, 'cov_rad': 1.67, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 1, 'num': 55, 'rgb': (0.05, 0.05, 0.5)}, 
'Ba': {'bond_ord_rad': 1.98, 'cov_rad': 1.34, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 2, 'num': 56, 'rgb': (0.05, 0.05, 0.5)}, 
'La': {'bond_ord_rad': 1.69, 'cov_rad': 1.87, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 12, 'num': 57, 'rgb': (0.05, 0.05, 0.5)}, 
'Ce': {'bond_ord_rad': 1.83, 'cov_rad': 1.83, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 58, 'rgb': (0.05, 0.05, 0.5)}, 
'Pr': {'bond_ord_rad': 1.82, 'cov_rad': 1.82, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 59, 'rgb': (0.05, 0.05, 0.5)}, 
'Nd': {'bond_ord_rad': 1.81, 'cov_rad': 1.81, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 60, 'rgb': (0.05, 0.05, 0.5)}, 
'Pm': {'bond_ord_rad': 1.8, 'cov_rad': 1.8, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 61, 'rgb': (0.05, 0.05, 0.5)}, 
'Sm': {'bond_ord_rad': 1.8, 'cov_rad': 1.8, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 62, 'rgb': (0.05, 0.05, 0.5)}, 
'Eu': {'bond_ord_rad': 1.99, 'cov_rad': 1.99, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 63, 'rgb': (0.05, 0.05, 0.5)}, 
'Gd': {'bond_ord_rad': 1.79, 'cov_rad': 1.79, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 64, 'rgb': (0.05, 0.05, 0.5)}, 
'Tb': {'bond_ord_rad': 1.76, 'cov_rad': 1.76, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 65, 'rgb': (0.05, 0.05, 0.5)}, 
'Dy': {'bond_ord_rad': 1.75, 'cov_rad': 1.75, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 66, 'rgb': (0.05, 0.05, 0.5)}, 
'Ho': {'bond_ord_rad': 1.74, 'cov_rad': 1.74, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 67, 'rgb': (0.05, 0.05, 0.5)}, 
'Er': {'bond_ord_rad': 1.73, 'cov_rad': 1.73, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 68, 'rgb': (0.05, 0.05, 0.5)}, 
'Tm': {'bond_ord_rad': 1.72, 'cov_rad': 1.72, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 69, 'rgb': (0.05, 0.05, 0.5)}, 
'Yb': {'bond_ord_rad': 1.94, 'cov_rad': 1.94, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 70, 'rgb': (0.05, 0.05, 0.5)}, 
'Lu': {'bond_ord_rad': 1.72, 'cov_rad': 1.72, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 71, 'rgb': (0.05, 0.05, 0.5)}, 
'Hf': {'bond_ord_rad': 1.44, 'cov_rad': 1.57, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 72, 'rgb': (0.05, 0.05, 0.5)}, 
'Ta': {'bond_ord_rad': 1.34, 'cov_rad': 1.43, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 73, 'rgb': (0.05, 0.05, 0.5)}, 
'W': {'bond_ord_rad': 1.3, 'cov_rad': 1.37, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 74, 'rgb': (0.05, 0.05, 0.5)}, 
'Re': {'bond_ord_rad': 1.28, 'cov_rad': 1.35, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 75, 'rgb': (0.05, 0.05, 0.5)}, 
'Os': {'bond_ord_rad': 1.26, 'cov_rad': 1.37, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 76, 'rgb': (0.05, 0.05, 0.5)}, 
'Ir': {'bond_ord_rad': 1.27, 'cov_rad': 1.32, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 77, 'rgb': (0.05, 0.05, 0.5)}, 
'Pt': {'bond_ord_rad': 1.3, 'cov_rad': 1.5, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 78, 'rgb': (0.05, 0.05, 0.5)}, 
'Au': {'bond_ord_rad': 1.34, 'cov_rad': 1.5, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 79, 'rgb': (0.05, 0.05, 0.5)}, 
'Hg': {'bond_ord_rad': 1.49, 'cov_rad': 1.7, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 80, 'rgb': (0.05, 0.05, 0.5)}, 
'Tl': {'bond_ord_rad': 1.48, 'cov_rad': 1.55, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 3, 'num': 81, 'rgb': (0.05, 0.05, 0.5)}, 
'Pb': {'bond_ord_rad': 1.48, 'cov_rad': 1.54, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 4, 'num': 82, 'rgb': (0.05, 0.05, 0.5)}, 
'Bi': {'bond_ord_rad': 1.45, 'cov_rad': 1.54, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 3, 'num': 83, 'rgb': (0.05, 0.05, 0.5)}, 
'Po': {'bond_ord_rad': 1.46, 'cov_rad': 1.68, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 2, 'num': 84, 'rgb': (0.05, 0.05, 0.5)}, 
'At': {'bond_ord_rad': 1.45, 'cov_rad': 0.7, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 1, 'num': 85, 'rgb': (0.05, 0.05, 0.5)}, 
'Rn': {'bond_ord_rad': 2.4, 'cov_rad': 2.4, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 0, 'num': 86, 'rgb': (0.05, 0.05, 0.5)}, 
'Fr': {'bond_ord_rad': 2.0, 'cov_rad': 2.0, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 1, 'num': 87, 'rgb': (0.05, 0.05, 0.5)}, 
'Ra': {'bond_ord_rad': 1.9, 'cov_rad': 1.9, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 2, 'num': 88, 'rgb': (0.05, 0.05, 0.5)}, 
'Ac': {'bond_ord_rad': 1.88, 'cov_rad': 1.88, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 89, 'rgb': (0.05, 0.05, 0.5)}, 
'Th': {'bond_ord_rad': 1.79, 'cov_rad': 1.79, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 90, 'rgb': (0.05, 0.05, 0.5)}, 
'Pa': {'bond_ord_rad': 1.61, 'cov_rad': 1.61, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 91, 'rgb': (0.05, 0.05, 0.5)}, 
'U': {'bond_ord_rad': 1.58, 'cov_rad': 1.58, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 92, 'rgb': (0.05, 0.05, 0.5)}, 
'Np': {'bond_ord_rad': 1.55, 'cov_rad': 1.55, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 93, 'rgb': (0.05, 0.05, 0.5)}, 
'Pu': {'bond_ord_rad': 1.53, 'cov_rad': 1.53, 'bs_rad': 0.18, 'vdw_rad': 0.9, 'max_bonds': 6, 'num': 94, 'rgb': (0.05, 0.05, 0.5)}, 
'Am': {'bond_ord_rad': 1.07, 'cov_rad': 1.51, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 95, 'rgb': (0.05, 0.05, 0.5)}, 
'Cm': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 96, 'rgb': (0.05, 0.05, 0.5)}, 
'Bk': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 97, 'rgb': (0.05, 0.05, 0.5)}, 
'Cf': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 98, 'rgb': (0.05, 0.05, 0.5)}, 
'Es': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 99, 'rgb': (0.05, 0.05, 0.5)}, 
'Fm': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 100, 'rgb': (0.05, 0.05, 0.5)}, 
'Md': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 101, 'rgb': (0.05, 0.05, 0.5)}, 
'No': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 102, 'rgb': (0.05, 0.05, 0.5)}, 
'Lr': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 106, 'rgb': (0.05, 0.05, 0.5)}, 
'D': {'bond_ord_rad': 0.346, 'cov_rad': 0.23, 'bs_rad': 0.24, 'vdw_rad': 1.2, 'max_bonds': 1, 'num': 107, 'rgb': (1.0, 1.0, 1.0)}, 
}

#
# guess_babel_elements
#
# These atom types are based on the assumption that the atom name in the molecule
# is only 1 character, and these assignments are the most likely equivalent
# if the second letter is missing.  This dictionary should be consulted only if 
# the previous dictionary babel_elements failed to contain the desired element.
#
guess_babel_elements = {
'X': {'bond_ord_rad': 0.0, 'cov_rad': 0.0, 'bs_rad': 0.0, 'vdw_rad': 0.0, 'max_bonds': 0, 'num': 0, 'rgb': (0.07, 0.5, 0.7)}, 
'L': {'bond_ord_rad': 1.23, 'cov_rad': 0.68, 'bs_rad': 0.304, 'vdw_rad': 1.52, 'max_bonds': 1, 'num': 3, 'rgb': (0.0, 0.0, 0.1)}, 
'T': {'bond_ord_rad': 1.32, 'cov_rad': 1.47, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 22, 'rgb': (0.05, 0.05, 0.5)}, 
'V': {'bond_ord_rad': 1.22, 'cov_rad': 1.33, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 23, 'rgb': (0.05, 0.05, 0.5)}, 
'M': {'bond_ord_rad': 1.17, 'cov_rad': 1.35, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 8, 'num': 25, 'rgb': (0.05, 0.05, 0.5)}, 
'F': {'bond_ord_rad': 1.17, 'cov_rad': 1.34, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 26, 'rgb': (0.7, 0.0, 0.1)}, 
'Z': {'bond_ord_rad': 1.25, 'cov_rad': 1.45, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 30, 'rgb': (0.05, 0.05, 0.5)}, 
'G': {'bond_ord_rad': 1.26, 'cov_rad': 1.22, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 3, 'num': 31, 'rgb': (0.05, 0.05, 0.5)}, 
'R': {'bond_ord_rad': 2.16, 'cov_rad': 1.47, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 1, 'num': 37, 'rgb': (0.05, 0.05, 0.5)}, 
'E': {'bond_ord_rad': 1.99, 'cov_rad': 1.99, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 63, 'rgb': (0.05, 0.05, 0.5)}, 
'A': {'bond_ord_rad': 1.34, 'cov_rad': 1.5, 'bs_rad': 0.34, 'vdw_rad': 1.7, 'max_bonds': 6, 'num': 79, 'rgb': (0.05, 0.05, 0.5)}, 
'D': {'bond_ord_rad': 0.346, 'cov_rad': 0.23, 'bs_rad': 0.24, 'vdw_rad': 1.2, 'max_bonds': 1, 'num': 107, 'rgb': (1.0, 1.0, 1.0)}, 
}
#END
