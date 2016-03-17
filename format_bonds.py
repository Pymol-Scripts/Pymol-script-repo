'''
format_bonds.py
Described at: http://www.pymolwiki.org/format_bonds
Version 1.0 (2014)
##################################################
Formats bonds in aromatic or charged residues
##################################################
Plugin contributed by Andreas Warnecke
(andreas.warnecke@ki.se, 4ndreas.warneck3@gmail.com)
##################################################
VERSION NOTES:
    1.0    2014    First release
'''

################################################################################
from __future__ import print_function
from pymol import cmd
from pymol import stored


def format_bonds(
    selection='all',
    bonds=4,
):
    '''
DESCRIPTION
    Formats bonds in aromatic or charged residues
EXAMPLE
    frag PHE
    format_bonds
USAGE
    format_bonds [ selection [, bonds ]]
ARGUMENTS
    selection: <str> input selection {default: 'all'}
    bonds:     <int> toogles format of bonds
               1: single bonds (deactivates valence display)
               2: regular double bonds (activates valence display)
             >=3: delocalized (activates valence display)
    '''
    # Selection
    try:
        # group selection with bracketing and select complete residues
        selection = '(byres (' + str(selection) + '))'
        # checks functional selection
        cmd.count_atoms(selection)
    except:
        print("invalid selection")
        return False

    # PARAMETERS
    try:
        bonds = int(bonds)
    except:
        pass
    if (not (bonds in [1, 2])):
        bonds = 4

    if bonds == 1:
        cmd.set('valence', 0)
        print("Valence display disabled!")
        return bonds
    else:
        cmd.set('valence', 1)
        print("Valence display enabled!")
    # proceed

    ##### SELECTION BY OBJECT AND CHAIN #####
    # variable for the selections
    # get the names of the proteins in the selection
    objects = cmd.get_object_list(selection)
    # include chains

    # subselect chains
    names = []
    for p in objects:
        for chain in cmd.get_chains('model ' + p) or ['']:
            names.append("(model %s and chain '%s')" % (p, chain))

    ##### SELECTION LISTS #####
    # get TRP
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn TRP+NIW) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    # the integer is to ensure unique keys
    TRP_tuple = (1,) + tuple(stored.temp)

    # get PHETYR
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn PHE+TYR+PTR+NIY+PNIY) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    PHETYR_tuple = (2,) + tuple(stored.temp)

    # get HIS
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn HIS) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    HIS_tuple = (3,) + tuple(stored.temp)

    # get NITRO
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn NIY+PNIY+NIW) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    NITRO_tuple = (4,) + tuple(stored.temp)

    # get GLU
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn GLU) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    GLU_tuple = (5,) + tuple(stored.temp)

    # get ASP
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn ASP) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    ASP_tuple = (6,) + tuple(stored.temp)

    # get CTERM
    stored.temp = []
    for p in names:
        cmd.iterate(
            '(byres (last %s)) and (not (hetatm)) '
            'and (name OXT)' % (p),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    CTERM_tuple = (7,) + tuple(stored.temp)

    # get ARG
    stored.temp = []
    for p in names:
        cmd.iterate((
            '(%s) and (resn ARG) '
            'and (name CA)' % (p)),
            'stored.temp.append("(%s and resi "+str(resi)+")")' % p
        )
    ARG_tuple = (8,) + tuple(stored.temp)

    ##### SELECTION TUPLES DONE #####

    ##### ATOM LISTS #####

    TRP_bonds_all = [
        ['CG', 'CD1'],
        ['CD1', 'NE1'],
        ['NE1', 'CE2'],
        ['CE2', 'CD2'],
        ['CD2', 'CG'],
        ['CD2', 'CE3'],
        ['CE3', 'CZ3'],
        ['CZ3', 'CH2'],
        ['CH2', 'CZ2'],
        ['CZ2', 'CE2']
    ]

    TRP_bonds_double = [
        ['CG', 'CD1'],
        ['CE2', 'CD2'],
        ['CE3', 'CZ3'],
        ['CH2', 'CZ2']
    ]

    PHETYR_bonds_all = [
        ['CG', 'CD1'],
        ['CD1', 'CE1'],
        ['CE1', 'CZ'],
        ['CZ', 'CE2'],
        ['CE2', 'CD2'],
        ['CD2', 'CG']
    ]

    PHETYR_bonds_double = [
        ['CG', 'CD1'],
        ['CE1', 'CZ'],
        ['CE2', 'CD2']
    ]

    HIS_bonds_all = [
        ['CG', 'CD2'],
        ['CD2', 'NE2'],
        ['NE2', 'CE1'],
        ['CE1', 'ND1'],
        ['ND1', 'CG'],
    ]

    HIS_bonds_double = [
        ['CG', 'CD2'],
        ['CE1', 'ND1']
    ]

    NITRO_bonds_all = [
        ['NN', 'O1'],
        ['NN', 'O2']
    ]
    NITRO_bonds_double = [
        ['NN', 'O1']
    ]

    GLU_bonds_all = [
        ['CD', 'OE1'],
        ['CD', 'OE2']
    ]
    GLU_bonds_double = [
        ['CD', 'OE1']
    ]

    ASP_bonds_all = [
        ['CG', 'OD1'],
        ['CG', 'OD2']
    ]
    ASP_bonds_double = [
        ['CG', 'OD1']
    ]

    CTERM_bonds_all = [
        ['C', 'O'],
        ['C', 'OXT']
    ]
    CTERM_bonds_double = [
        ['C', 'O']
    ]

    ARG_bonds_all = [
        ['CZ', 'NH1'],
        ['CZ', 'NH2']
    ]

    ARG_bonds_double=[
    ['CZ','NH1']
    ]

    ##### FORMATING #####

    # dictionary: entries:atoms
    format_dict = {
        TRP_tuple: [TRP_bonds_all, TRP_bonds_double],
        PHETYR_tuple: [PHETYR_bonds_all, PHETYR_bonds_double],
        HIS_tuple: [HIS_bonds_all, HIS_bonds_double],
        NITRO_tuple: [NITRO_bonds_all, NITRO_bonds_double],
        GLU_tuple: [GLU_bonds_all, GLU_bonds_double],
        ASP_tuple: [ASP_bonds_all, ASP_bonds_double],
        CTERM_tuple: [CTERM_bonds_all, CTERM_bonds_double],
        ARG_tuple: [ARG_bonds_all, ARG_bonds_double]
    }

    if bonds != 2:
        lines = 4
        print("Formating as delocalized bonds")
    else:
        lines = 1
        print("Formating as double bonds")

    # for all tuples (i.e format_dict.keys())
    for p in list(format_dict.keys()):
        # go through list except ID at pos 1
        for q in p[1:]:
            # format bonds
            for r in format_dict[p][0]:
                cmd.unbond('%s and name %s' % (q, r[0]), '%s and name %s' % (q, r[1]))
                cmd.bond('%s and name %s' % (q, r[0]), '%s and name %s' % (q, r[1]), lines)
            if lines == 1:
                # add double bonds
                for r in format_dict[p][1]:
                    cmd.unbond('%s and name %s' % (q, r[0]), '%s and name %s' % (q, r[1]))
                    cmd.bond('%s and name %s' % (q, r[0]), '%s and name %s' % (q, r[1]), 2)

    return bonds
cmd.extend("format_bonds", format_bonds)
cmd.auto_arg[0]['format_bonds'] = [lambda: cmd.Shortcut(['all', 'resn PHE+TYR+PTR+NIY+PNIY+TRP+NIW+HIS', 'resn GLU+ASP', 'resn ARG', 'last all']), 'selection=', ', ']
cmd.auto_arg[1]['format_bonds'] = [lambda: cmd.Shortcut(['4', '2', '1']), 'bonds=', '']
################################################################################
