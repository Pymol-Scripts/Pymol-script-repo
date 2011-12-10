#!/usr/bin/env python
"""
Calulate all torsion angles and assign secondary structure and mesostate code.

This module uses much of the code from the original BIOMOL collection of
utilities written by Raj Srinivasani with enhancements by Nick Fitzkee.

The script was put together by Pat Fleming so that a user would not need
to have the BIOMOL distribution installed to run PROSS.

Note that since Raj's time the definitions of mesostates has been superceded
by the fine grained 30 deg x 30 deg grid for most purposes. Either mesostate
grid will work for PROSS. Give your choice as an argument (see USAGE below).

Date: September 2004
Author: Pat Fleming, pat.fleming@jhu.edu
"""

USAGE = """
python PROSS.py pdbfile [oldmeso | fgmeso]
Choose old mesostate definitions or finegrained mesostate definitions.
Default is finegrained.
"""

import sys
import string, re
import copy
import gzip
import math
import types

# This script does not require Numeric
HAVE_NUMPY = 0    

_RAD_TO_DEG = 180.0/math.pi
_DEG_TO_RAD = math.pi/180.0

RESIDUES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN', 'GLU',
            'GLX', 'GLY', 'HIS', 'HIP', 'ILE', 'LEU', 'LYS', 'MET',
            'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ACE',
            'NME', 'PCA', 'FOR', 'ASX', 'AMD']
ONE_LETTER_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'LYS': 'K'}

_CHI1_DEFAULT = ('N', 'CA', 'CB', 'CG')
CHI1_ATOMS = {
    'ARG': _CHI1_DEFAULT,
    'ASN': _CHI1_DEFAULT,
    'ASP': _CHI1_DEFAULT,
    'CYS': ('N', 'CA', 'CB', 'SG'),
    'SER': ('N', 'CA', 'CB', 'OG'),
    'GLN': _CHI1_DEFAULT,
    'GLU': _CHI1_DEFAULT,
    'HIS': _CHI1_DEFAULT,
    'ILE': ('N', 'CA', 'CB', 'CG1'),
    'LEU': _CHI1_DEFAULT,
    'LYS': _CHI1_DEFAULT,
    'MET': _CHI1_DEFAULT,
    'PHE': _CHI1_DEFAULT,
    'PRO': _CHI1_DEFAULT,
    'THR': ('N', 'CA', 'CB', 'OG1'),
    'TRP': _CHI1_DEFAULT,
    'TYR': _CHI1_DEFAULT,
    'VAL': ('N', 'CA', 'CB', 'CG1'),
    'GLX': _CHI1_DEFAULT,
    'ASX': _CHI1_DEFAULT,
    'CYX': _CHI1_DEFAULT,
    'PCA': _CHI1_DEFAULT
    }

CHI2_ATOMS = {
    'ARG': ('CA', 'CB', 'CG', 'CD'),
    'ASN': ('CA', 'CB', 'CG', 'OD1'),
    'ASP': ('CA', 'CB', 'CG', 'OD1'),
    'ASX': ('CA', 'CB', 'CG', 'OD1'),
    'GLN': ('CA', 'CB', 'CG', 'CD'),
    'GLU': ('CA', 'CB', 'CG', 'CD'),
    'GLX': ('CA', 'CB', 'CG', 'CD'),
    'HIS': ('CA', 'CB', 'CG', 'ND1'),
    'ILE': ('CA', 'CB', 'CG1', 'CD1'),
    'LEU': ('CA', 'CB', 'CG', 'CD1'),
    'LYS': ('CA', 'CB', 'CG', 'CD'),
    'MET': ('CA', 'CB', 'CG', 'SD'),
    'PHE': ('CA', 'CB', 'CG', 'CD1'),
    'PRO': ('CA', 'CB', 'CG', 'CD'),
    'TRP': ('CA', 'CB', 'CG', 'CD1'),
    'TYR': ('CA', 'CB', 'CG', 'CD1')
    }

CHI3_ATOMS = {
    'ARG': ('CB', 'CG', 'CD', 'NE'),
    'GLN': ('CB', 'CG', 'CD', 'OE1'),
    'GLU': ('CB', 'CG', 'CD', 'OE1'),
    'GLX': ('CB', 'CG', 'CD', 'OE1'),
    'LYS': ('CB', 'CG', 'CD', 'CE'),
    'MET': ('CB', 'CG', 'SD', 'CE')
    }

CHI4_ATOMS = {
    'ARG': ('CG', 'CD', 'NE', 'CZ'),
    'LYS': ('CG', 'CD', 'CE', 'NZ')
    }


# Unpack PDB Line Constants - used for simple index changes
UP_SERIAL     = 0
UP_NAME       = 1
UP_ALTLOC     = 2
UP_RESNAME    = 3
UP_CHAINID    = 4
UP_RESSEQ     = 5
UP_ICODE      = 6
UP_X          = 7
UP_Y          = 8
UP_Z          = 9
UP_OCCUPANCY  = 10
UP_TEMPFACTOR = 11

##########################
# Start Mesostate Definition Code 
##########################

# This is redundant but kept for compatibility
default = 'fgmeso'

MSDEFS = {}

# Setup Fine Grain Mesostate Bins (ala Pat Fleming)

MSDEFS['fgmeso'] = {}
mydefs = MSDEFS['fgmeso']

mydefs['RC_DICT'] = {( 180,-165): 'Aa', ( 180,-135): 'Ab', ( 180,-105): 'Ac',
                     ( 180, -75): 'Ad', ( 180, -45): 'Ae', ( 180, -15): 'Af',
                     ( 180,  15): 'Ag', ( 180,  45): 'Ah', ( 180,  75): 'Ai',
                     ( 180, 105): 'Aj', ( 180, 135): 'Ak', ( 180, 165): 'Al',
                     (-150,-165): 'Ba', (-150,-135): 'Bb', (-150,-105): 'Bc',
                     (-150, -75): 'Bd', (-150, -45): 'Be', (-150, -15): 'Bf',
                     (-150,  15): 'Bg', (-150,  45): 'Bh', (-150,  75): 'Bi',
                     (-150, 105): 'Bj', (-150, 135): 'Bk', (-150, 165): 'Bl',
                     (-120,-165): 'Ca', (-120,-135): 'Cb', (-120,-105): 'Cc',
                     (-120, -75): 'Cd', (-120, -45): 'Ce', (-120, -15): 'Cf',
                     (-120,  15): 'Cg', (-120,  45): 'Ch', (-120,  75): 'Ci',
                     (-120, 105): 'Cj', (-120, 135): 'Ck', (-120, 165): 'Cl',
                     ( -90,-165): 'Da', ( -90,-135): 'Db', ( -90,-105): 'Dc',
                     ( -90, -75): 'Dd', ( -90, -45): 'De', ( -90, -15): 'Df',
                     ( -90,  15): 'Dg', ( -90,  45): 'Dh', ( -90,  75): 'Di',
                     ( -90, 105): 'Dj', ( -90, 135): 'Dk', ( -90, 165): 'Dl',
                     ( -60,-165): 'Ea', ( -60,-135): 'Eb', ( -60,-105): 'Ec',
                     ( -60, -75): 'Ed', ( -60, -45): 'Ee', ( -60, -15): 'Ef',
                     ( -60,  15): 'Eg', ( -60,  45): 'Eh', ( -60,  75): 'Ei',
                     ( -60, 105): 'Ej', ( -60, 135): 'Ek', ( -60, 165): 'El',
                     ( -30,-165): 'Fa', ( -30,-135): 'Fb', ( -30,-105): 'Fc',
                     ( -30, -75): 'Fd', ( -30, -45): 'Fe', ( -30, -15): 'Ff',
                     ( -30,  15): 'Fg', ( -30,  45): 'Fh', ( -30,  75): 'Fi',
                     ( -30, 105): 'Fj', ( -30, 135): 'Fk', ( -30, 165): 'Fl',
                     (   0,-165): 'Ja', (   0,-135): 'Jb', (   0,-105): 'Jc',
                     (   0, -75): 'Jd', (   0, -45): 'Je', (   0, -15): 'Jf',
                     (   0,  15): 'Jg', (   0,  45): 'Jh', (   0,  75): 'Ji',
                     (   0, 105): 'Jj', (   0, 135): 'Jk', (   0, 165): 'Jl',
                     (  30,-165): 'Ha', (  30,-135): 'Hb', (  30,-105): 'Hc',
                     (  30, -75): 'Hd', (  30, -45): 'He', (  30, -15): 'Hf',
                     (  30,  15): 'Hg', (  30,  45): 'Hh', (  30,  75): 'Hi',
                     (  30, 105): 'Hj', (  30, 135): 'Hk', (  30, 165): 'Hl',
                     (  60,-165): 'Ia', (  60,-135): 'Ib', (  60,-105): 'Ic',
                     (  60, -75): 'Id', (  60, -45): 'Ie', (  60, -15): 'If',
                     (  60,  15): 'Ig', (  60,  45): 'Ih', (  60,  75): 'Ii',
                     (  60, 105): 'Ij', (  60, 135): 'Ik', (  60, 165): 'Il',
                     (  90,-165): 'Ja', (  90,-135): 'Jb', (  90,-105): 'Jc',
                     (  90, -75): 'Jd', (  90, -45): 'Je', (  90, -15): 'Jf',
                     (  90,  15): 'Jg', (  90,  45): 'Jh', (  90,  75): 'Ji',
                     (  90, 105): 'Jj', (  90, 135): 'Jk', (  90, 165): 'Jl',
                     ( 120,-165): 'Ka', ( 120,-135): 'Kb', ( 120,-105): 'Kc',
                     ( 120, -75): 'Kd', ( 120, -45): 'Ke', ( 120, -15): 'Kf',
                     ( 120,  15): 'Kg', ( 120,  45): 'Kh', ( 120,  75): 'Ki',
                     ( 120, 105): 'Kj', ( 120, 135): 'Kk', ( 120, 165): 'Kl',
                     ( 150,-165): 'La', ( 150,-135): 'Lb', ( 150,-105): 'Lc',
                     ( 150, -75): 'Ld', ( 150, -45): 'Le', ( 150, -15): 'Lf',
                     ( 150,  15): 'Lg', ( 150,  45): 'Lh', ( 150,  75): 'Li',
                     ( 150, 105): 'Lj', ( 150, 135): 'Lk', ( 150, 165): 'Ll'}

# What mesostate to use when the angles are invalid (e.g. 999.99)
mydefs['INVALID'] = '??'

# What mesostate to use when the omega angle is odd (e.g. < 90.0)
mydefs['OMEGA'] = '**'

# The size of mesostate codes used in this set.
mydefs['CODE_LENGTH'] = 2

# Geometric parameters: DELTA is the size of the bins, XXX_OFF is
# the offset from 0 degrees of the centers of the bins.

mydefs['DELTA']   = 30.0 
mydefs['PHI_OFF'] =  0.0
mydefs['PSI_OFF'] = 15.0

# Set up turns and turn dictionary.  Dictionary contains the number of
# occurences in the PDB.  This number isn't used, so for new turn types
# '1' is sufficient (as is used for PII, below)

mydefs['TURNS'] = {'EfDf': 5226, 'EeEf': 4593, 'EfEf': 4061, 'EfDg': 3883,
                   'EeDg': 2118, 'EeEe': 1950, 'EfCg': 1932, 'EeDf': 1785,
                   'EkJf': 1577, 'EkIg': 1106, 'EfEe':  995, 'EkJg':  760,
                   'EeCg':  553, 'DfDf':  479, 'EfCf':  395, 'DgDf':  332,
                   'DfDg':  330, 'IhIg':  310, 'EfDe':  309, 'EkIh':  298,
                   'DgCg':  275, 'DfCg':  267, 'IbDg':  266, 'DfEe':  260,
                   'FeEf':  250, 'IbEf':  249, 'DfEf':  219, 'IhJf':  216,
                   'IhJg':  213, 'IgIg':  207, 'EfCh':  188, 'DgEe':  180,
                   'DgEf':  176, 'EeEg':  172, 'IhIh':  153, 'EeDe':  150,
                   'IgJg':  147, 'EkKf':  147, 'EeCh':  147, 'IbDf':  131,
                   'DgDg':  128, 'EgDf':  127, 'FeDg':  114, 'ElIg':  111,
                   'IgIh':  107, 'DfDe':  107, 'EjIg':  101, 'EeCf':  100,
                   'DfCh':   94, 'DgCf':   91, 'DfCf':   91, 'DeEe':   91,
                   'DkIh':   88, 'FeDf':   79, 'EkIf':   78, 'EeDh':   76,
                   'DgCh':   74, 'IgJf':   71, 'EjJg':   71, 'FeEe':   69,
                   'DlIh':   66, 'EgCg':   65, 'ElIh':   62, 'EjJf':   62,
                   'FeCg':   59, 'DlIg':   56, 'IbCg':   54, 'EfEg':   54,
                   'EkJe':   53, 'FkJf':   52, 'ElJg':   51, 'DgDe':   49,
                   'DlJg':   46, 'EgCf':   45, 'IaEf':   40, 'FkIg':   39,
                   'JaEf':   38, 'EjIh':   38, 'EgEf':   38, 'DkJg':   36,
                   'DeEf':   34, 'EeCi':   31, 'JgIh':   29, 'IcEf':   29,
                   'EkKe':   29, 'DkIg':   29, 'IbEe':   27, 'EgDg':   27,
                   'EeFe':   27, 'EjKf':   26, 'IaDf':   25, 'HhIg':   24,
                   'HbDg':   24, 'ElJf':   24, 'EfDh':   24, 'IcDf':   23,
                   'EfBh':   23, 'IcDg':   22, 'IcCg':   22, 'FkJg':   21,
                   'FeCh':   21, 'IgKf':   20, 'FdDg':   20, 'EkHh':   20,
                   'DfDh':   20, 'DgBh':   19, 'DfBh':   19, 'DeDf':   19,
                   'DfFe':   18, 'EfFe':   17, 'EgEe':   16, 'EgDe':   16,
                   'DkJf':   16, 'JgJg':   15, 'IbEg':   15, 'IbCh':   15,
                   'EfBg':   15, 'DgCe':   15, 'JlEf':   14, 'CgCg':   14,
                   'HhJf':   13, 'EeBi':   13, 'DfBi':   13, 'IhIf':   12,
                   'FeEg':   12, 'FdEf':   12, 'EdEf':   12, 'DlJf':   12,
                   'DhCg':   12, 'JgIg':   11, 'IeBg':   11, 'FjIg':   11,
                   'FdCh':   11, 'EdEe':   11, 'JfIh':   10, 'JaEe':   10,
                   'HhJg':   10, 'HbEf':   10, 'HbCh':   10, 'FkIh':   10,
                   'FjJf':   10, 'ElJe':   10, 'DhDf':   10, 'CgDf':   10}

# Set up the PII defitions, similar to dictionary above.

mydefs['PII'] = {'Dk':1, 'Dl':1, 'Ek':1, 'El':1}

# Set up the codes that define helix and strand.  Here, rather than storing
# the dictionary like we did above, we'll store the regular expression
# matcher directly.  This prevents us from recompiling it every time we
# want to find a helix or strand.

helix  = ('De', 'Df', 'Ed', 'Ee', 'Ef', 'Fd', 'Fe')
strand = ('Bj', 'Bk', 'Bl', 'Cj', 'Ck', 'Cl', 'Dj', 'Dk', 'Dl')

pat_helix  = "(%s){5,}" % string.join(map(lambda x: "(%s)" % x, helix), '|')
pat_strand = "(%s){3,}" % string.join(map(lambda x: "(%s)" % x, strand), '|')

mydefs['HELIX'] = re.compile(pat_helix)
mydefs['STRAND'] = re.compile(pat_strand)

###########################
# Setup Old Mesostate Bins (ala Raj Srinivasan)
##########################

MSDEFS['oldmeso'] = {}
mydefs = MSDEFS['oldmeso']

mydefs['RC_DICT'] = {( 180, 180): 'A', ( 180,-120): 'B', ( 180, -60): 'C',
                     ( 180,   0): 'D', ( 180,  60): 'E', ( 180, 120): 'F',
                     (-120, 180): 'G', (-120,-120): 'H', (-120, -60): 'I',
                     (-120,   0): 'J', (-120,  60): 'K', (-120, 120): 'L',
                     ( -60, 180): 'M', ( -60,-120): 'N', ( -60, -60): 'O',
                     ( -60,   0): 'P', ( -60,  60): 'Q', ( -60, 120): 'R',
                     (   0, 180): 'S', (   0,-120): 'T', (   0, -60): 'U',
                     (   0,   0): 'V', (   0,  60): 'W', (   0, 120): 'X',
                     (  60, 180): 'm', (  60,-120): 'r', (  60, -60): 'q',
                     (  60,   0): 'p', (  60,  60): 'o', (  60, 120): 'n',
                     ( 120, 180): 'g', ( 120,-120): 'l', ( 120, -60): 'k',
                     ( 120,   0): 'j', ( 120,  60): 'i', ( 120, 120): 'h'}

# What mesostate to use when the angles are invalid (e.g. 999.99)
mydefs['INVALID'] = '?'

# What mesostate to use when the omega angle is odd (e.g. < 90.0)
mydefs['OMEGA'] = '*'

# The size of mesostate codes used in this set.
mydefs['CODE_LENGTH'] = 1

# Geometric parameters: DELTA is the size of the bins, XXX_OFF is
# the offset from 0 degrees of the centers of the bins.

mydefs['DELTA']   = 60.0 
mydefs['PHI_OFF'] =  0.0
mydefs['PSI_OFF'] =  0.0

# Set up turns and turn dictionary.  Dictionary contains the type
# of the turn (no primes)

mydefs['TURNS'] = {'OO': 1, 'OP': 1, 'OJ': 1, 'PO': 1, 'PP': 1, 'PJ': 1,
                   'JO': 1, 'JP': 1, 'JJ': 1,
                   'Mo': 2, 'Mp': 2, 'Mj': 2, 'Ro': 2, 'Rp': 2, 'Rj': 2,
                   'oo': 3, 'op': 3, 'oj': 3, 'po': 3, 'pp': 3, 'pj': 3,
                   'jo': 3, 'jp': 3, 'jj': 3,
                   'mO': 4, 'mP': 4, 'mJ': 4, 'rO': 4, 'rP': 4, 'rJ': 4}

# Set up the PII defitions, similar to dictionary above.

mydefs['PII'] = {'M':1, 'R':1}

# Set up the codes that define helix and strand.  Here, rather than storing
# the dictionary like we did above, we'll store the regular expression
# matcher directly.  This prevents us from recompiling it every time we
# want to find a helix or strand.

helix  = ('O', 'P')
strand = ('L', 'G', 'F', 'A')

pat_helix  = "(%s){5,}" % string.join(map(lambda x: "(%s)" % x, helix),  '|')
pat_strand = "(%s){3,}" % string.join(map(lambda x: "(%s)" % x, strand), '|')

mydefs['HELIX'] = re.compile(pat_helix)
mydefs['STRAND'] = re.compile(pat_strand)

##########################
# End Mesostate Definition Code 
##########################

def res_rc(r1, r2, r3=180, mcodes=None):
    """res_rc(r1, r2, r3) - get mesostate code for a residue

    Given a phi (r1), psi (r2), and omega (r3) torsion angle, calculate
    the mesostate code that describes that residue. 
    A mesostate will be returned in
    all but two cases:  if omega deviates from planarity by more than 90
    degrees, '*' is returned.  Also, if any torsions are greater than
    180.0 (biomol assignes 999.0 degrees to angles that that are
    indeterminate), prosst.INVALID is returned.  Here, r3 defaults to
    180.0.
    """

    if not mcodes: mcodes = default
    ms = MSDEFS[mcodes]
    OMEGA   = ms['OMEGA']
    INVALID = ms['INVALID']
    PHI_OFF = ms['PHI_OFF']
    PSI_OFF = ms['PSI_OFF']
    DELTA   = ms['DELTA']
    RC_DICT = ms['RC_DICT']

    if (abs(r3) <= 90.0):
        return OMEGA
    elif r1>180.0 or r2>180.0 or r3>180.0:
        return INVALID

    ir1 = -int(PHI_OFF) + int(round((r1+PHI_OFF)/DELTA )) * int(DELTA)
    ir2 = -int(PSI_OFF) + int(round((r2+PSI_OFF)/DELTA )) * int(DELTA)

    while ir1 <= -180: ir1 = ir1 + 360
    while ir1 >   180: ir1 = ir1 - 360
    while ir2 <= -180: ir2 = ir2 + 360
    while ir2 >   180: ir2 = ir2 - 360

    return RC_DICT[(ir1,ir2)]

def rc_codes(chain, phi=None, psi=None, ome=None, mcodes=None):
    """rc_codes(chain, phi, psi, ome) - return rotamer codes

    Given a protein chain (and optionally phi, psi, omega), this
    function will return a list of mesostate codes that
    applies to the chain, as determined by res_rc.
    """
    if not mcodes: mcodes = default
    n = range(len(chain))
    if phi is None: phi = map(chain.phi, n)
    if psi is None: psi = map(chain.psi, n)
    if ome is None: ome = map(chain.omega, n)
    return map(lambda x, y, z: res_rc(x, y, z, mcodes), phi, psi, ome)

def rc_ss(chain, phi=None, psi=None, ome=None, mcodes=None):
    """rc_ss(chain, phi, psi, ome) - calculate secondary structure

    This function calculates the secondary structure using the PROSS method
    with rotamer codes.  Given a chain, and optionally a list of phi,
    psi, and omega, it calculates the backbone secondary structure of
    the chain.  The return value is (phi, psi, ome, sst), where
    phi, psi, and ome are calculated if not specified, and sst is the
    secondary structure codes: H = helix, E = strand, P = PII, C = coil.
    """

    if not mcodes: mcodes = default
    ms = MSDEFS[mcodes]
    PII    = ms['PII']
    TURNS  = ms['TURNS']
    HELIX  = ms['HELIX']
    STRAND = ms['STRAND']

    if phi is None:
        chain.gaps()
        nres = len(chain)
        phi = map(chain.phi, xrange(nres))
    else:
        nres = len(chain)
    if psi is None: psi = map(chain.psi, xrange(nres))
    if ome is None: ome = map(chain.omega, xrange(nres))

    codes = rc_codes(chain, phi, psi, ome, mcodes)

    chain.gaps()

    sst = ['C']*nres

    is_PII = PII.has_key

    for i in xrange(nres-1):
        code = codes[i]
        if is_PII(code):
            sst[i] = 'P'

    is_turn = TURNS.has_key

    for i in xrange(nres-1):
        code = codes[i] + codes[i+1]
        if is_turn(code):
            sst[i] = sst[i+1] = 'T'

    helices = _rc_find(codes, HELIX, mcodes)
    strands = _rc_find(codes, STRAND, mcodes)

    for helix in helices:
        i, j = helix

        for k in range(i, j):
            sst[k] = 'H'

    for strand in strands:
        i, j = strand
        for k in range(i, j):
            if sst[k] in ('C', 'P'): sst[k] = 'E'
#            if sst[k] == 'C': sst[k] = 'E'

    return phi, psi, ome, sst
    
def _rc_find(codes, pattern, mcodes=None):
    """_rc_find(codes, pat_obj) - find a endpoints of a regexp

    Given a list of mesostate codes, this function identifies a endpoints
    of a match  <pattern>.  <pat_obj> is a compiled regular expression
    pattern whose matches will be returned as pairs indicated start,
    end in <codes>
    """
    if not mcodes: mcodes = default
    CODE_LENGTH = MSDEFS[mcodes]['CODE_LENGTH']

    if not type(codes) == type(''):
        codes = string.join(codes, '')

    matches = []
    it = pattern.finditer(codes)

    try:
        while 1:
            mat = it.next()
            matches.append((mat.start()/CODE_LENGTH, mat.end()/CODE_LENGTH))
    except StopIteration:
        pass

    return matches

            
##############################
# Protein objects and methods
##############################

def is_atom(o):
    return hasattr(o, '_is_an_atom3d')


def is_residue(o):
    return hasattr(o, '_is_a_residue')

def is_chain(o):
    return hasattr(o, '_is_a_chain')

def is_mol(o):
    return hasattr(o, '_is_a_mol')

####################


HAVE_POP = hasattr([], 'pop')
HAVE_EXTEND = hasattr([], 'extend')

class TypedList:
    """A Python list restricted to having objects of the same type.
    An instance of a TypedList is created as follows:
    
    mylist = TypedList(function, [elements])
    
    where function is a python function which takes an argument
    and returns 1 or 0 indicating whether the object represented
    by the argument is of the correct type, and elements is an optional
    list of elements to be added into the instance.  Here is a
    full blown example.

    def is_int(o):
        return type(o) == type(0)

    mylist = TypedList(is_int, [0, 1, 2, 3])

    New elements are added to the list as follows:
    mylist.append(25)

    Instances of TypedList support all operations available for
    Python Lists (as of Python version 1.5.2a2)
    """

    _is_a_typed_list = 1
    
    def __init__(self, function, elements=None):
        self._func = function
        self.elements = []
        if not elements is None:
            if self._func(elements):
                self.elements.append(elements)
            else:
                for el in elements:
                    self.append(el)

    def append(self, el):
        if self._func(el):
            self.elements.append(el)
        else:
            raise TypeError, 'Element to be added to list has incorrect type.'

    def __len__(self): return len(self.elements)

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__, self._func.__name__,
                              `self.elements`)

    def __str__(self):
        return `self.elements`
    
    def __getitem__(self, i): return self.elements[i]
    
    def __setitem__(self, i, v):
        if self._func(v):
            self.elements[i] = v
        else:
            raise TypeError, 'Item not of correct type in __setitem__'

    def __delitem__(self, i): del self.elements[i]

    def __getslice__(self, i, j):
        new = self.__class__(self._func)
        new.elements = self.elements[i:j]
        return new

    def __setslice__(self, i, j, v):
        if self._alltrue(v):
            self.elements[i:j] = v

    def __delslice__(self, i, j):
        del self.elements[i:j]

    def __add__(self, other):
        if not hasattr(other, '_is_a_typed_list'):
            raise TypeError,'List to be concatenated not instance of %s' %\
                  self.__class__
        if self._func <> other._func:
            raise TypeError, 'Lists to be added not of same type'
        new = self.__class__(self._func)
        new.elements = self.elements + other.elements
        return new

    def __mul__(self, other):
        if type(other) == type(0):
            new = self.__class__(self._func)
            new.elements = self.elements * other
            return new
        else:
            raise TypeError, "can't multiply list with non-int"
        
    __rmul__ = __mul__

    def __copy__(self):
        new = self.__class__(self._func)
        for el in self.elements:
            new.elements.append(el.__copy__())
        have = new.__dict__.has_key
        for key in self.__dict__.keys():
            if not have(key):
                new.__dict__[key] = copy.deepcopy(self.__dict__[key])
        return new

    __deepcopy__ = clone = __copy__

    def _alltrue(self, els):
        return len(els) == len(filter(None, map(self._func, els)))

    def sort(self): self.elements.sort()

    def reverse(self): self.elements.reverse()
    
    def count(self, el): return self.elements.count(el)
    
    def extend(self, els):
        if self._alltrue(els):
            if HAVE_EXTEND:
                self.elements.extend(els)
            else:
                for el in els:
                    self.elements.append(el)
        else:
            raise TypeError, 'One or more elements of list not of correct type'
                    
    def pop(self):
        if HAVE_POP:
            return self.elements.pop()
        else:
            el = self.elements[-1]
            del self.elements[-1]
            return el
    
    def index(self, el): return self.elements.index(el)
    
    def remove(self, el): self.elements.remove(el)
    
    def insert(self, pos, el):
        if self._func(el):
            self.elements.insert(pos, el)
        else:
            raise TypeError, 'Item not of correct type in insert'

    def indices(self, len=len):
        return xrange(len(self.elements))

    def reverse_indices(self, len=len):
        return xrange(len(self.elements)-1, -1, -1)

    
    
        
        

####################

class molResidue(TypedList):
    _is_a_residue = 1

    def __init__(self, name='', atoms=None, **kw):
        TypedList.__init__(self, is_atom, atoms)
        self.name = name
        for key, value in kw.items():
            setattr(self, key, value)        
        
    def num_atoms(self):
        """returns the number of atoms in residue"""
        return len(self.elements)

    def has_atom(self, name):
        """returns true if residue has an atom named 'name'"""
        for atom in self.elements:
            if atom.name == name:
                return 1
        return 0

    def atoms(self):
        return self.elements

    def atoms_with_name(self, *names):
        """returns atoms in residue with specified names"""
        ret = []
        Append = ret.append
        for name in names:
            for atom in self.elements:
                if atom.name == name:
                    Append(atom)
                    break
        return ret

    def delete_atoms_with_name(self, *names):
        """delete atoms in residue with specified names"""
        els = self.elements
        for i in self.reverse_indices():
            atom = els[i]
            if atom.name in names:
                del els[i]

    def atoms_not_with_name(self, *names):
        """returns atoms in residue excluding specified names"""
        ret = []
        for atom in self.elements:
            if not atom.name in names:
                ret.append(atom)
        return ret

    def atom_coordinates(self, *names):
        """returns coordinates of named atoms.
        If names is omitted all atom coordinates are returned."""

        if len(names)==0:
            atoms = self.elements
        else:
            atoms = apply(self.atoms_with_name, names)

        na = len(atoms)
        if na == 0: return

        if HAVE_NUMPY:
            a = Numeric.zeros((na, 3), 'd')
        else:
            a = [None]*na

        pos = 0
        for atom in atoms:
            a[pos] = atom.coords()
            pos = pos + 1
        return a

    def assign_radii(self):
        raise AttributeError, 'Should be defined in specialized class'

    def type(self):
        return 'residue'

class molChain(TypedList):
    _is_a_chain = 1

    def __init__(self, name='', residues=None, **kw):
        self.name = name
        TypedList.__init__(self, is_residue, residues)
        for key, value in kw.items():
            setattr(self, key, value)        


    def num_residues(self): return len(self)

    def num_atoms(self):
        na = 0
        for res in self.elements:
            na = na + len(res.elements)
        return na

    def atoms(self):
        ret = []
        Append = ret.append
        for res in self.elements:
            for atom in res.elements:
                Append(atom)
        return ret

    def residues_with_name(self, *names):
        """returns named residues as a python list"""
        if len(names) == 0:
            return
        l = []
        for res in self.elements:
            if res.name in names:
                l.append(res)
        return l

    def delete_residues_with_name(self, *names):
        """delete named residues from Chain"""
        if len(names) == 0:
            return
        els = self.elements
        for i in self.reverse_indices():
            if els[i].name in names:
                del els[i]

    def residues_not_with_name(self, *names):
        """returns residues excluding specified names as a python list"""
        ret = []
        for res in self.elements:
            if not res.name in names:
                ret.append(res)
        return ret
        

    def atoms_with_name(self, *names):
        ret = []
        Append = ret.append
        if len(names) > 0:
            for res in self.elements:
                for name in names:
                    for atom in res.elements:
                        if atom.name == name:
                            Append(atom)
                            break
        else:
            for res in self.elements:
                for atom in res.elements:
                    Append(atom)
                    
        return ret
            
    def delete_atoms_with_name(self, *names):
        """delete atoms in residue with specified names"""
        for res in self.elements:
            apply(res.delete_atoms_with_name, names)

    def atoms_not_with_name(self, *names):
        """returns atoms in residue excluding specified names"""
        ret = []
        for res in self.elements:
            ret[len(ret):] = apply(res.atoms_not_with_name, names)
        return ret

    def atom_coordinates(self, *names):
        """returns coordinates of named atoms. if names is None
        all atom coordinates are returned."""

        coords = []
        if len(names) > 0:
            for res in self.elements:
                for atom in res.elements:
                    if atom.name in names:
                        coords.append(atom.coords())
        else:
            atoms = apply(self.atoms_with_name, names)
            coords = map(lambda a:a.coords(), atoms)
        if HAVE_NUMPY:
            return Numeric.array(coords)
        else:
            return coords

    def atoms(self):
        atoms = []
        for res in self.elements:
            for atom in res: atoms.append(atom)
        return atoms

    def delete_alt_locs(self):
        """delete_alt_locs - remove secondary conformations in the chain

    In a chain with multiple occupancy and alternate location identifiers,
    it is often desirable to eliminate the secondary conformations for
    use in simulation, etc.  This function (abitrarily) finds and selects
    the first given conformation and deletes all other conformations.
        """
        AtomCount = self.present
        chain = self.elements
        
        delete = []

        for i in xrange(len(chain)):
            residue = chain[i]
            rid = (residue.idx, residue.icode)

            for j in xrange(len(residue)):
                atom = residue[j]
                anam = atom.name

                try:
                    acnt = AtomCount[rid][anam]
                except KeyError:
                    print "Unable to locate %s %s %s in present dictionary."%\
                          (rid[0], rid[1], anam)
                    return

                if acnt == 1:
                    continue
                if acnt < 1:
                    AtomCount[rid][anam] = acnt + 1
                    delete.append((i, j))
                    continue

                atom.alt = ' '
                AtomCount[rid][anam] = -acnt + 2

        delete.reverse()
        for r, a in delete:
            del chain[r][a]

    def assign_radii(self):
        for res in self:
            res.assign_radii()

    def preserve_chain_hetero(self):
        """prevent hetero residues from being deleted as hetero atoms

    Normally, delete_hetero will delete all hetero atoms from a
    molecule. This includes waters and heterogroups (hemes, etc.),
    but it also includes hetero residues--nonstandard residues
    that have backbone connectivity but perhaps extra atoms (e.g.
    S-hydroxy-cysteine).  Deleting these residues may disrupt an
    otherwise continuous chain and may be undesirable.
    
    Given that a chain has a valid SEQRES entry, this function will
    'unset' the hetero flag for heterogroups that are involved in
    the sequence itself.  When delete_hetero is run, these atoms
    will be preserved.
    """
        for i in xrange(self.num_residues()):
            if hasattr(self[i], 'chain_het') and hasattr(self[i], 'het'):
                delattr(self[i], 'het')
                
    def delete_hetero(self):
        """delete all hetero atoms from a protein

    This function removes all HETATM records from the molecular
    structure.  This include waters, heme groups, as well as residues
    with nonstandard structures
    """
        for i in xrange(self.num_residues()-1, -1, -1):
            if hasattr(self[i], 'het'):
                del self[i]

    def delete_waters(self):
        for i in xrange(self.num_residues()-1, -1, -1):
            if self[i].name == 'HOH':
                del self[i]

    def translate(self, dx, dy=None, dz=None):
        for res in self.elements:
            res.translate(dx, dy, dz)
    
    def rotatex(self, theta):
        for res in self.elements:
            res.rotatex(theta)

    def rotatey(self, theta):

            res.rotatey(theta)

    def rotatez(self, theta):
        for res in self.elements:
            res.rotatez(theta)

    def type(self):
        return 'Chain'

class molMol(TypedList):
    _is_a_mol = 1

    def __init__(self, name='', chains=None):
        self.name = name
        self.resolution  = None
        self.method = []
        self.rprog = None
        self.rfree = None
        self.rvalu = None
        TypedList.__init__(self, is_chain, chains)

    def num_chain(self): return len(self)

    def num_res(self):
        nr = 0
        for chain in self: nr = nr + len(chain)
        return nr

    def num_atoms(self):
        na = 0
        for chain in self: na = na + chain.num_atoms()
        return na

    def atoms(self):
        ret = []
        Append = ret.append

        for chain in self.elements:
            for res in chain.elements:
                for atom in res.elements:
                    Append(atom)
        return ret

    def chains_with_name(self, *names):
        chains = []
        for chain in self.elements:
            if chain.name in names: chains.append(chain)
        return chains

    def residues_with_name(self, *names):
        residues = []
        for chain in self.elements:
            for res in chain.elements:
                if res.name in names: residues.append(res)
        return residues


    def atoms_with_name(self, *names):
        ret = []
        Append = ret.append
        for chain in self.elements:
            for res in chain.elements:
                for name in names:
                    for atom in res.elements:
                        if atom.name == name:
                            Append(atom)
                            break
        return ret

    def delete_atoms_with_name(self, *names):
        """delete atoms in assembly with specified names"""
        for chain in self.elements:
            apply(chain.delete_atoms_with_name, names)

    def atoms_not_with_name(self, *names):
        """returns atoms in assembly excluding specified names"""
        ret = []
        for chain in self.elements:
            ret[len(ret):] = apply(chain.atoms_not_with_name, names)
        return ret


    def atom_coordinates(self, *names):

        coords = []
        if len(names) == 0:
            for chain in self.elements:
                for res in chain.elements:
                    for atom in res.elements:
                        coords.append(atom.coords())
        else:
            atoms = apply(self.atoms_with_name, names)
            coords = map(lambda a:a.coords(), atoms)
            del atoms
        if HAVE_NUMPY:
            return Numeric.array(coords)
        else:
            return coords

    def delete_alt_locs(self):
        """delete_alt_locs - remove all secondary conformations"""
        for chain in self:
            chain.delete_alt_locs()

    def assign_radii(self):
        for chain in self:
            chain.assign_radii()

    def delete_water(self):
        for chain in self:
            chain.delete_water()

    def preserve_chain_hetero(self):
        for i in self.indices():
            self[i].preserve_chain_hetero()

    def delete_hetero(self):
        for i in self.reverse_indices():
            self[i].delete_hetero()
            if len(self[i]) == 0:
                del self[i]

    def pdb(self, file, renum=0, seq=1):
        """pdb(file, renum=0, seq=1) - write a PDB file

    Given a file name/open file handle (with .gz option), this function
    will write the PDB coordinate file corresponding to the current
    molecular object.  If renum is true, then residue indices will be
    renumbered starting from 1, if seq is true, then SEQRES records
    will be written in the PDB file.

    This function is a wrapper for pdbout.write_pdb.
    """
        pdbout.write_pdb(self, file, renum, seq)
                            


def pack_pdb_line(atom, idx, aan, aanum, ic, cn, f, fmt):
    """pack_pdb_line - construct an ATOM record for a PDB file

    pack_pdb_line takes the following arguments:

    o  atom   An atom3d object, used for extracting coordinates, occupancy
              atom name, alternate location, and B-factor.  If these
              attributes aren't contained in the object, reasonable defaults
              are used.

    o  idx    The atom serial number

    o  aan    The name of the containing amino residue (DNA base)

    o  aanum  The number of the containing residue (DNA base)

    o  ic     The insertion code for homologous residue numbering

    o  cn     The name of the chain

    o  f      An open file handle where the format will be written.

    o  fmt    A format to be used.  Generally this is either for ATOM or
              HETATM records.
    """

    try:
        atn = pad(atom.name)
    except AttributeError:
        atn = " X  "
    
    try:
        al = atom.alt
    except AttributeError:
        al = ' '

    try:
        occ = atom.occ
        bf  = atom.bf
    except AttributeError:
        occ = 1.00
        bf  = 25.00

    x = atom.xcoord()
    y = atom.ycoord()
    z = atom.zcoord()

#    fmt = '%s %s %s %s %s %s %s %s %s %s %s\n'

#    atfmt = 'ATOM  %5i %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n'
#    htfmt = 'HETATM%5i %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n'

    f.write(fmt % (idx, atn, al, aan, cn, aanum, ic, x, y, z, occ, bf))

def pad(nm):
    """pad(nm) - a function to pad an atom name with appropraiate spaces"""
    space = '    '

    if len(nm) >= 4:
        return nm

    try:
        int(nm[0])
    except ValueError:
        nm = ' ' + nm

    return nm + space[0:4-len(nm)]


def write_pdb(myMol, f, renum=0, seq=None, pack=pack_pdb_line):
    """write_pdb - write a PDB file from a molecular object

    write_pdb takes a biomol.Mol object and creates a PDB file
    based on the definitions given.  It can support all atomic
    properties as well as insertion codes and alternate locations.
    Essentially, biomol is designed such that if the original molecule
    has all that stuff it will be retained in the final PDB file.  If
    that is not desired, you must make the modifications to the Mol
    object yourself (or use one of the utilities).

    write_pdb takes the following arguments:

    o  myMol   The **mol.Mol** object to be written.

    o  f       A file instruction.  This can be a string, in which
               case f will be interpreted as a filename.  If the
               string f ends with .gz, the resulting file will be
               compressed using gzip.  Alternatively, f can be an
               open file object (it must have the write attribute).
               If f is an object, the PDB will be written to the
               object.  When f is a file object, SEQRES and END
               records will not be written by default.

    o  seq     A boolean that determines whether SEQRES records
               will be written.  By default this is true (except in
               the case where f is an object).
    """

    if not is_mol(myMol):
        print "This function only works for biomol.mol types."
        return

    default_seq = 1
    
    if hasattr(f, 'write'):
        out = f
        close = 0
    elif f[-3:] == '.gz':
        import gzip
        out = gzip.open(f, 'w')
        close = 1
    else:
        out = open(f, 'w')
        default_seq = 0
        close = 1

    if seq is None:
        seq = default_seq

    if seq:
        for chain in myMol:
            write_sequences(chain, out)

    atfmt = 'ATOM  %5i %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n'
    htfmt = 'HETATM%5i %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n'

    snum = 1
    for chain in myMol:
        cname = chain.name
        rname = ''
        rnum = 0
        rlast = int(chain[0].idx)-1
        
        for res in chain:
            rname = res.name
            icode = res.icode
            fmt = atfmt
            if hasattr(res, 'het') or hasattr(res, 'chain_het'):
                fmt = htfmt
                
            if not renum:
                rnum = int(res.idx)
            else:
                rnum = rnum + (int(res.idx) - rlast)
                icode = ' '

            for atom in res:
                pack(atom, snum, rname, rnum, icode, cname, out, fmt)
                snum = snum + 1

            if hasattr(res, 'gap'):
                terminate(snum, rname, cname, rnum, icode, out)
                snum = snum + 1

            rlast = int(res.idx)

        terminate(snum, rname, cname, rnum, icode, out)
        snum = snum + 1

    if close:
        out.write('END\n')
        out.close()

        

def write_sequences(chain, f):
    """write_sequences - write SEQRES records to a PDB file

    write_sequences takes the following arguments:

    o  chain   A mol.Chain type whose residues/bases will be displayed

    o  f       An open file handle where information will be written.
    """
    
    nr = len(chain)
    nm = chain.name
    sn = 1
    elem = 0
    fmt = 'SEQRES  %2i %1s %4i  ' 

    for i in xrange(nr):
        if elem == 0:
            f.write(fmt % (sn, nm, nr))
            sn = sn + 1

        f.write('%3s ' % chain[i].name)
        elem = elem + 1

        if elem == 13:
            f.write('\n')
            elem = 0

    if elem != 0:
        f.write('\n')


def terminate(snum, rn, cn, rnum, icode, f):
    """terminate - write a TER record to a PDB file"""
    fmt = 'TER   %5i      %3s %1s%4s%1s\n'
    f.write(fmt % (snum, rn, cn, rnum, icode))


class Atom3d:
    _is_an_atom3d = 1

    def __init__(self, x, y, z, cnf={}):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        if cnf:
            for key, value in cnf.items():
                self.__dict__[key] = value

    def __repr__(self):
        return "Atom3d(%s, %s, %s)" % (self.x, self.y, self.z)
    __str__ = __repr__

    def coords(self):
        if HAVE_NUMPY:
            return Numeric.array((self.x, self.y, self.z))
        else:
            return [self.x, self.y, self.z]

    def set_coords(self, x, y=None, z=None):
        if y is None:
            self.x, self.y, self.z = map(float, x)
        else:
            self.x, self.y, self.z = float(x), float(y), float(z)

    def xcoord(self, v=None):
        if v is None:
            return self.x
        else:
            self.x = float(v)

    def ycoord(self, v=None):
        if v is None:
            return self.y
        else:
            self.y = float(v)

    def zcoord(self, v=None):
        if v is None:
            return self.z
        else:
            self.z = float(v)

    def clone(self):
        new = Atom3d(0., 0., 0.)
        new.__dict__.update(self.__dict__)
        return new

    __copy__ = clone

#   def __deepcopy__(self):
#       cnf = copy.deepcopy(self.__dict__)
#       return self.__class__(0, 0, 0, cnf)

    def distance(self, other, sqrt=math.sqrt):
        if not hasattr(other, '_is_an_atom3d'):
            raise TypeError, 'distance argument must be Atom3d type'
        
        return sqrt( (self.x - other.x)**2 +
                     (self.y - other.y)**2 +
                     (self.z - other.z)**2 )

    def sqr_distance(self, other):
        if not hasattr(other, '_is_an_atom3d'):
            raise TypeError, 'sqr_distance argument must be Atom3d type'

        return (self.x - other.x)**2 + \
               (self.y - other.y)**2 + \
               (self.z - other.z)**2

    def angle(self, a, b, sqrt=math.sqrt, acos=math.acos):
        if not hasattr(a, '_is_an_atom3d') or \
           not hasattr(b, '_is_an_atom3d'):
            raise TypeError, 'angle arguments must be Atom3d type'

        x1, y1, z1 = self.x - a.x, self.y - a.y, self.z - a.z
        x2, y2, z2 = b.x - a.x, b.y - a.y, b.z - a.z
        v11 = x1**2 + y1**2 + z1**2
        v22 = x2**2 + y2**2 + z2**2
        
        if (v11 == 0.0 or v22 == 0.0):
            raise ValueError, 'Null vector in angle'
        
        v12 = x1*x2 + y1*y2 + z1*z2
        
        ang = v12/sqrt(v11*v22)
        
        if ang >= 1.0:
            return 0.0
        elif ang <= -1.0:
            return 180.0
        else:
            return acos(ang) * _RAD_TO_DEG

    def torsion(self, a, b, c, sqrt=math.sqrt, acos=math.acos):
        
        if not hasattr(a, '_is_an_atom3d') or \
           not hasattr(b, '_is_an_atom3d') or \
           not hasattr(c, '_is_an_atom3d'):
            raise TypeError, 'torsion arguments must be Atom3d type'
        
        v12x, v12y, v12z = self.x - a.x, self.y - a.y, self.z - a.z
        v32x, v32y, v32z = b.x - a.x, b.y - a.y, b.z - a.z
        v43x, v43y, v43z = c.x - b.x, c.y - b.y, c.z - b.z

        vn13x = v12y*v32z - v12z*v32y
        vn13y = v12z*v32x - v12x*v32z
        vn13z = v12x*v32y - v12y*v32x

        vn24x = v32z*v43y - v32y*v43z
        vn24y = v32x*v43z - v32z*v43x
        vn24z = v32y*v43x - v32x*v43y

        v12 = vn13x*vn24x + vn13y*vn24y + vn13z*vn24z
        v11 = vn13x**2 + vn13y**2 + vn13z**2
        v22 = vn24x**2 + vn24y**2 + vn24z**2

        ang = v12/sqrt(v11*v22)
        if ang >= 1.0:
            return 0.0
        elif ang <= -1.0:
            return -180.0
        else:
            ang = acos(ang) * _RAD_TO_DEG
            
        vtmp = vn13x * (vn24y*v32z - vn24z*v32y) + \
               vn13y * (vn24z*v32x - vn24x*v32z) + \
               vn13z * (vn24x*v32y - vn24y*v32x) < 0.0
        if vtmp:
            return -ang
        else:
            return ang

########################
# functions that operate on a collection of atoms 
########################


def fromint(v1, dis, v2, a, v3, t, sqrt=math.sqrt, sin=math.sin, cos=math.cos):

    ang = a * _DEG_TO_RAD
    sina = sin(ang)
    cosa = cos(ang)

    tors = t * _DEG_TO_RAD 
    sint = sina * sin(tors)
    cost = sina * cos(tors)

    x1, y1, z1 = v1.coords()
    x2, y2, z2 = v2.coords()
    x3, y3, z3 = v3.coords()
    
    u1x = x2 - x3
    u1y = y2 - y3
    u1z = z2 - z3
    d = 1.0 / sqrt(u1x*u1x + u1y*u1y + u1z*u1z)
    u1x = u1x*d
    u1y = u1y*d
    u1z = u1z*d

    u2x = x1 - x2
    u2y = y1 - y2
    u2z = z1 - z2;
    d = 1.0 / sqrt(u2x*u2x + u2y*u2y + u2z*u2z)
    u2x = u2x*d
    u2y = u2y*d
    u2z = u2z*d

    cosine = u1x*u2x + u1y*u2y + u1z*u2z

    if (abs(cosine) < 1.0):
        sine = 1.0/sqrt(1.0 - cosine*cosine)
    else:
        sine = 1.0/sqrt(cosine*cosine - 1.0)
            
    u3x = sine * (u1y*u2z - u1z*u2y)
    u3y = sine * (u1z*u2x - u1x*u2z)
    u3z = sine * (u1x*u2y - u1y*u2x)

    u4x = cost * (u3y*u2z - u3z*u2y)
    u4y = cost * (u3z*u2x - u3x*u2z)
    u4z = cost * (u3x*u2y - u3y*u2x)

    return Atom3d(x1 + dis*(-u2x*cosa + u4x + u3x*sint),
                  y1 + dis*(-u2y*cosa + u4y + u3y*sint),
                  z1 + dis*(-u2z*cosa + u4z + u3z*sint))

def _atof(s, atof=string.atof):
    try:
        return atof(s)
    except:
        return None

def _atoi(s, atoi=string.atoi):
    try:
        return atoi(s)
    except:
        return None

def get_sequences(file):
    if file[-3:] == '.gz':
        ff = gzip.GzipFile(file)
    else:
        ff = open(file, 'r')
    
    # read until we find a line with SEQRES card
    line = ff.readline()
    while line:
        if line[:6] == 'SEQRES':
            break
        else:
            line = ff.readline()

    if not line: # no SEQRES records
        return {}

    sequences = {}
    chain = line[11]
    sequences[chain] = []
    
    while line[:6] == 'SEQRES':
        chain = line[11]
        residues = string.split(line[19:71])
        if not sequences.has_key(chain):
            sequences[chain] = []
        sequences[chain].extend(residues)
        line = ff.readline()

    return sequences
     

def unpack_pdb_line(line, ATOF=_atof, ATOI=_atoi, STRIP=string.strip):
    return (ATOI(line[6:11]),
           STRIP(line[12:16]),
           line[16],
           STRIP(line[17:21]),
           line[21],
           STRIP(line[22:26]),
           line[26],                     # Insertion of residues (?)
           ATOF(line[30:38]),
           ATOF(line[38:46]),
           ATOF(line[46:54]),
           ATOF(line[54:60]),
           ATOF(line[60:66]))


def atom_build(t, atom=Atom3d):
    atm = atom(t[UP_X], t[UP_Y], t[UP_Z]);
    atm.occ = t[UP_OCCUPANCY];
    atm.bf = t[UP_TEMPFACTOR];
    atm.name = t[UP_NAME]; 
    atm.idx = t[UP_SERIAL];
    atm.alt = t[UP_ALTLOC];
    return atm

def _type_mol(mol):
    if not len(mol):
        return mol
    for i in range(len(mol)):
        chain = mol[i]
        if is_protein(chain):
            Chain = Protein
            Residue = AminoAcid
        else:
            continue

        new_chain = Chain(chain.name, present=chain.present, model=chain.model)
        Append = new_chain.elements.append
        for res in chain.elements:
            new_res = Residue(res.name, idx=res.idx, icode=res.icode)
            new_res.elements = res.elements[:]
            Append(new_res)
            if hasattr(res, 'het'):
                new_res.het = 1
            if hasattr(res, 'chain_het'):
                new_res.chain_het = 1
        mol[i] = new_chain
        del chain
    return mol

def is_protein(chain):
    for res in chain.elements:
        if res.name in RESIDUES:
            return 1
        #elif hasattr(res, 'het') and not hasattr(res, 'chain_het'):
        #    return 0
    return 0

###########################
# End of protein objects and  functions
###########################

def read_pdb(f, as_protein=0, as_rna=0, as_dna=0, all_models=0,
             unpack=unpack_pdb_line, atom_build=atom_build):
    """read_pdb - read the entire contents of a PDB file as a molcule object

    This function parses a PDB file and returns a heirarchic data structure
    that can be used to browse atomic contents.  The heirarchic data elements
    are built on 'typed' lists.  The first level is a list of 'chains',
    where a chain can be non-specific or a protein, DNA, or RNA chain.  A
    chain (of any type) is a list of residues, and a residue is a list of
    atoms.  Atoms are defined as a full-blown Python object, with its own
    member functions and data (including x, y, z).  Residues and chains are
    derived from the TypedList object, but can have their own member data
    and functions--for example, the *pross* function is only available for
    Protein chains, and chains in general have a public variable *name*
    the described the SegID for that chain.

    read_pdb takes the following arguments:

    o  f           A file containing PDB data.  If *f* is a string, it will
                   be checked for a '.gz' extension.  Given that it contains
                   GZipped data, the file will be opened as a GZipFile and
                   parsed accordingly.  If no '.gz' extension is present,
                   the file will be opened as a regular text file.
                   Alternatively, f can be an open file object with the
                   'readline' member function.  If this is the case, the file
                   will not be closed upon exiting, and read_pdb will read
                   until the 'END' token is reached.

    o  as_protein  Read the PDB as a protein.  It will have the protein
                   specific member functions that allow calculation of torsion
                   angles, etc.  By default this is false, and the type is
                   determined automatically from the residue types.

    o  as_rna      Read the PDB as RNA.  As above, this is false by default
                   and determined automatically.

    o  as_dna      Read the PDB as DNA.  As above, this is false by default
                   and determined automatically.

    o  all_models  A boolean.  The default is false.  If true, all NMR models
                   will be read in as differing chains rather than just the
                   first one.
    """
    
    if type(f) is type(''):
        close = 1
        if f[-3:] == '.gz':
            inp = gzip.GzipFile(f)
        else:
            inp = open(f, 'r')
    else:
        inp = f
        close = 0

    try:
        sequences = get_sequences(inp.name)
    except AttributeError:
        sequences = get_sequences(inp.filename)

    Mol = molMol
    Chain = molChain
    Residue = molResidue

    if as_protein:
        Chain = Protein
        Residue = AminoAcid

    new_mol = Mol('')
    
    nextline = inp.readline
    while 1:
        line = nextline()
        if not line:
            break

        # Try to determine the resolution, in angstroms.
        try:
            if len(line) > 10 and line[:6] == 'REMARK' and line[9] == '2':
                line = line[10:70].upper()
                rpos = string.find(line, 'RESOLUTION.')
                if rpos >= 0:
                    new_mol.resolution = float(line[rpos+11:].split()[0])
        except:
            pass

        # Try to pick out expdata information, but don't attempt to
        # parse the specific type.
        try:
            if line[:6] == 'EXPDTA':
                for meth in line[10:].split(';'):
                    new_mol.method.append(meth.strip())
        except:
            pass

        # Try to identify refinement parameters: Program name,
        # free-R, and standard refinement factor.
        try:
            if len(line) > 10 and line[:6] == 'REMARK' and line[9] == '3':
                line = line[10:70].upper()

                if string.find(line, ':') >= 0:
                    get_val = lambda x: string.split(x, ':')[1].strip()
                else:
                    get_val = lambda x: string.split(x)[0]
                    
                ppos = string.find(line,  '  PROGRAM  ')
                rfpos = string.find(line, '  FREE R VALUE  ')
                rvpos = string.find(line, '  R VALUE  ')

                if ppos >= 0:
                    new_mol.rprog = get_val(line[ppos+11:])
                if rfpos >= 0:
                    new_mol.rfree = float(get_val(line[rfpos+16:]))
                if rvpos >= 0:
                    new_mol.rvalu = float(get_val(line[rvpos+11:]))
        except:
            pass
                
        if line[:6] == 'COMPND':
            if string.find(line, 'MOL_ID') <> -1: continue
            new_mol.name = new_mol.name + string.strip(line[10:72]) + ' '
        if line[:4] == 'ATOM' or line[:4] == 'HETA' or line[:5] == 'MODEL':
            break
    new_mol.name = string.strip(new_mol.name)    
    data = inp.readlines()
    if line[:4] in ('ATOM', 'HETA', 'MODE'):
        data.insert(0, line)

    old_chain = None
    old_res = None
    old_icode = None
    modnum = 0

    for line in data:
        key = line[:4]
        if key == 'ATOM' or key == 'HETA':
            t = unpack(line)

            if t[UP_CHAINID] <> old_chain:
                old_chain = t[UP_CHAINID]
                new_chain = Chain(old_chain, present={}, model=modnum)
                AppendResidue = new_chain.elements.append
                AtomsPresent = new_chain.present
                new_mol.append(new_chain)
                old_res = None
                old_icode = None
                
            if t[UP_RESSEQ] <> old_res or t[UP_ICODE] <> old_icode:
                old_res = t[UP_RESSEQ]
                old_icode = t[UP_ICODE]
                new_res = Residue(t[UP_RESNAME], idx=old_res, icode=old_icode)
                #if key == 'HETA' and not t[UP_RESNAME] in \
                #       sequences.get(t[UP_CHAINID], (' ',)):
                #    new_res.het = 1

                # All hetero groups are marked as hetero, but hetero
                # residues that are part of the actual chain
                # connectivity are marked as chain_het.  These atoms
                # can be 'unmarked as hetero' using the preserve_chain_hetero
                # function.
                if key == 'HETA':
                    if t[UP_RESNAME] in sequences.get(t[UP_CHAINID], (' ',)):
                        new_res.chain_het = 1
                    new_res.het = 1
                    
                AppendAtom = new_res.elements.append
                AppendResidue(new_res)
                if not AtomsPresent.has_key((old_res, old_icode)):
                    AtomsPresent[(old_res, old_icode)] = {}

            if AtomsPresent[(old_res, old_icode)].has_key(t[UP_NAME]):
                AtomsPresent[(old_res, old_icode)][t[UP_NAME]] = \
                     AtomsPresent[(old_res, old_icode)][t[UP_NAME]] + 1
            else:
                AtomsPresent[(old_res, old_icode)][t[UP_NAME]] = 1
            
            AppendAtom(atom_build(t))

        elif key == 'MODE':
            if len(new_mol) > 0 and not all_models:
                break
            #new_chain = Chain(string.split(line)[1])
            #AppendResidue = new_chain.elements.append
            #new_mol.append(new_chain)
            #old_chain = ' '

            modnum = int(string.split(line)[1])
            old_chain = None
            old_res = None
            old_icode = None

    del data
    if close:
        inp.close()

    

    if as_protein or as_dna or as_rna:
        return new_mol
    else:
        return _type_mol(new_mol)

def _read_chain(f, as_protein=0, as_rna=0, as_dna=0, unpack=unpack_pdb_line,
                atom_build = atom_build):

    if as_protein:
        Chain = Protein
        Residue = AminoAcid
    else:
        Chain = molChain
        Residue = molResidue

    nextline = f.readline

    new_chain = Chain('')
    AppendResidue = new_chain.elements.append
    old_chain = ' '
    old_res = None
    old_icode = None

    while 1:
        line = nextline()
        if not line:
            break
        key = line[:4]
        if key == 'ATOM' or key == 'HETA':
            break
        elif key == 'MODE':
            new_chain.name = string.split(line)[1]
            line = None

    if line:
        t = unpack(line)
        new_res = Residue(t[UP_RESNAME], idx=t[UP_RESSEQ], icode=t[UP_ICODE])
        AppendResidue(new_res)
        old_res = t[UP_RESSEQ]
        old_icode = t[UP_ICODE]
        AppendAtom = new_res.elements.append
        AppendAtom(atom_build(t))
        ResidueHasAtom = new_res.atoms_with_name

    while 1:
        line = nextline()
        if not line:
            break
        key = line[:4]
        if key == 'ATOM' or key == 'HETA':
            t = unpack(line)
            if t[UP_RESSEQ] <> old_res or t[UP_ICODE] <> old_icode:
                old_res = t[UP_RESSEQ]
                old_icode = t[UP_ICODE]
                new_res = Residue(old_res, idx=t[UP_RESSEQ], icode=t[UP_ICODE])
                AppendResidue(new_res)
                AppendAtom = new_res.elements.append
                ResidueHasAtom = new_res.atoms_with_name
            if not ResidueHasAtom(t[UP_NAME]):
                AppendAtom(atom_build(t))

        elif key[:3] in ('TER', 'END'):
            break

    if as_protein or as_dna or as_rna:
        return new_chain
    else:
        return _type_mol([new_chain])[0]
            
##########################


class Protein(molChain):
    _is_a_protein = 1

    def clean(self):
        """clean() - remove residues if they lack N, CA, or C atoms"""
        res = self.elements
        for i in self.reverse_indices():
            tmp = len(res[i].atoms_with_name('N', 'CA', 'C'))
            if tmp <> 3:
                del res[i]


    def gaps(self):
        """gaps() - identify gaps in the chain

    This function does a simple check to determine whether CA atoms are
    contiguous in the peptide chain.  First it, checks that CA atoms
    exist for all residues.  If that's okay, it checks to see if the
    distance is less than 18**0.5 A.  If either of these checks fail,
    residues involved in a gap have their .gap attribute set to true.
    """
        residues = self.elements
        dlist = []
        for i in xrange(len(self)-1):
            ca1 = residues[i].atoms_with_name('CA')
            if residues[i].name == 'ACE':
                ca1 = residues[i].atoms_with_name('CH3')
                
            ca2 = residues[i+1].atoms_with_name('CA')
            if residues[i+1].name == 'NME':
                ca2 = residues[i+1].atoms_with_name('CH3')
            
            if not ca1 or not ca2:
                if hasattr(residues[i+1], 'het'): continue
                self[i].gap = 1
            else:
                d = ca1[0].sqr_distance(ca2[0]) 
                if d > 18.0:
                    residues[i].gap = 1
                elif d <= 0.1:
                    dlist.append(i+1)
        dlist.reverse()
        for i in dlist:
            del residues[i]

    def phi(self, i):
        res = self.elements[i]       
        if i==0 or hasattr(self.elements[i-1], 'gap'):
            return 999.99
        
        try:
            a, b, c = res.atoms_with_name('N', 'CA', 'C')
        except:
            return 999.99

        try:
            if self[i-1].name == 'ACE':
                prev = self[i-1].atoms_with_name('CH3')[0]
            else:
                prev = self[i-1].atoms_with_name('C')[0]
        except IndexError:
            return 999.99
        else:
            return prev.torsion(a, b, c)

    def psi(self, i):
        res = self.elements[i]       
        if i >= len(self)-1 or hasattr(res, 'gap'):
            return 999.99

        try:
            a, b, c = res.atoms_with_name('N', 'CA', 'C')
        except:
            return 999.99

        try:
            next = self[i+1].atoms_with_name('N')[0]
        except IndexError:
            return 999.99
        else:
            return a.torsion(b, c, next)

    def omega(self, i):
        """omega(i) - calculate omega torsion for index i

    This function calcualtes the value of the omega torsion for
    index i.  Note: this value is widely calculated using atoms
    CA(i-1), C(i-1), N(i), and CA(i).  These atoms are used for
    this calculation.  However, other LINUS programs use the fllowing
    atoms: CA(i), C(i), N(i+1), CA(i+1).  Forewarned is forearmed.
    """
        res = self.elements[i]       

#        if i >= len(self)-1 or hasattr(res, 'gap'):
#            return 999.99

#        try:
#            a, b = res.atoms_with_name('CA', 'C')
#        except:
#            return 999.99

#        try:
#            c, d = self[i+1].atoms_with_name('N', 'CA')
#        except:
#            return 999.99

        if not i or hasattr(res, 'gap'):
            return 999.99

        try:
            if self[i-1].name == 'ACE':
                a, b = self[i-1].atoms_with_name('CH3', 'C')
            else:
                a, b = self[i-1].atoms_with_name('CA', 'C')
        except:
            return 999.99

        try:
            c, d = res.atoms_with_name('N', 'CA')
        except:
            return 999.99

        return a.torsion(b, c, d)

    def chi1(self, i):
        res = self.elements[i]       
        atoms = CHI1_ATOMS.get(res.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(res.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def chi2(self, i):
        res = self.elements[i]       
        atoms = CHI2_ATOMS.get(res.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(res.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def chi3(self, i):
        res = self.elements[i]       
        atoms = CHI3_ATOMS.get(res.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(res.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def chi4(self, i):
        res = self.elements[i]       
        atoms = CHI4_ATOMS.get(res.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(res.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def backbone_atoms(self, *ids):
        if not len(ids):
            ids = xrange(len(self))

        bb_atoms = ('N', 'CA', 'C', 'O')
        residues = []
        Append = residues.append
        res_class = res.__class__
        Strip = string.strip

        for i in ids:
            res = self[i]
            atoms = []
            for atom in res:
                if atom.name in bb_atoms: atoms.append(atom)
            new = res_class(res.name)
            res.elements = atoms
            Append(res)
        new_chain = self.__class__(self.name)
        new_chain.elements = residues
        return new_chain

    def torsions(self, i=None, map=map):
        if not i is None:
            return self.phi(i), self.psi(i), self.omega(i), self.chi1(i),\
                   self.chi2(i), self.chi3(i), self.chi4(i)
        else:
            n = range(len(self))
            return map(self.phi, n),\
                   map(self.psi, n),\
                   map(self.omega, n),\
                   map(self.chi1, n),\
                   map(self.chi2, n),\
                   map(self.chi3, n),\
                   map(self.chi4, n)


    def pross(self, phi=None, psi=None, ome=None, mcodes=None):
        return rc_ss(self, phi, psi, ome, mcodes)

    def codes(self, phi=None, psi=None, ome=None, mcodes=None):
        return rc_codes(self, phi, psi, ome, mcodes)
        
    def type(self):
        return 'protein'

    def sequence(self, one=0):
        n = self.num_residues()
        seq = map(getattr, self.elements, ('name',)*n)
        if one:
            get_one = ONE_LETTER_CODES.get
            seq = map(get_one, seq, 'X'*n)
        return seq
            

class AminoAcid(molResidue):

    def chi1(self):
        atoms = CHI1_ATOMS.get(self.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(self.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def chi2(self):
        atoms = CHI2_ATOMS.get(self.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(self.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def chi3(self):
        atoms = CHI3_ATOMS.get(self.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(self.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def chi4(self):
        atoms = CHI4_ATOMS.get(self.name, None)
        if atoms is None:
            return 999.99
        
        try:
            a, b, c, d = apply(self.atoms_with_name, atoms)
        except:
            return 999.99
        else:
            return a.torsion(b, c, d)

    def assign_radii(self):
        d_get = SASARAD[self.name].get
        for atom in self.elements:
            atom.radius = d_get(atom.name,0.0)
            

    
    

##########################



class PDBFile:

    def __init__(self, filename):
        self.filename = filename
        self._openfile()

    def _openfile(self):
        if self.filename[-3:] == '.gz':
            self.file = gzip.GzipFile(self.filename)
        else:
            self.file = open(self.filename)

    def read(self, as_protein=0, as_rna=0, as_dna=0):
        return read_pdb(self.file, as_protein, as_rna, as_dna)

    def read_chain(self, as_protein=0, as_rna=0, as_dna=0):
        return _read_chain(self.file, as_protein, as_rna, as_dna)



def run(file,mcode):
    mol = read_pdb(file)
    mol.delete_hetero()
    fmt = "%5s %3s  %s %s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n"

# Default is to screen, change here if file output wanted
#   angfile = open('angles.dat', 'w')
    angfile = sys.stdout

    for chain in mol:
        if chain.type() <> 'protein': continue
        chain.gaps()
        angfile.write('Chain: %s\n' % chain.name)
        angfile.write('          SS MS    phi      psi     omega     chi1     chi2     chi3     chi4\n')
        phi, psi, ome, chi1, chi2, chi3, chi4 = chain.torsions()
        phi, psi, ome, ss = chain.pross(phi, psi, ome)
        pos = 0
        for res in chain.elements:
            ms = res_rc(phi[pos], psi[pos],mcodes=mcode)
            angfile.write( fmt % (res.idx, res.name, ss[pos], ms, phi[pos], psi[pos], ome[pos],
                         chi1[pos], chi2[pos], chi3[pos], chi4[pos]))
            pos = pos + 1

#   angfile.close()

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 1:
        print USAGE
        sys.exit()

    # Default is finegrained mesostates
    if len(sys.argv) < 3:
        mcode = 'fgmeso'
    else:
        mcode = sys.argv[2]

    file = sys.argv[1]

    run(file, mcode)
