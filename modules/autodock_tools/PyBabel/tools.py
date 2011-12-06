#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/tools.py,v 1.1.1.1 2001/04/03 19:47:47 gillet Exp $
#
# $Id: tools.py,v 1.1.1.1 2001/04/03 19:47:47 gillet Exp $
#
#
#

import string

def read_element_table(filename):
    """void <- read_element_table(filename)
    populates the elementsTable dictionary from the a given file.
    the file provides:
    line number, element string, cov_rad, bond_ord_rad, vdw_rad, bs_rad,
    max_bonds, red, green, blue
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    elemTable = {}
    for i in range(len(lines)):
        dd = string.split(lines[i])
        elemTable[dd[1]] = { 'num':i,
                             'cov_rad':float(dd[2]),
                             'bond_ord_rad':float(dd[3]),
                             'vdw_rad':float(dd[4]),
                             'bs_rad':float(dd[5]),
                             'max_bonds':int(dd[6]),
                             'rgb': (float(dd[7]),float(dd[8]),float(dd[9]))
                           }
    return elemTable


def writeElementTableAsPythonCode(elemTab, inFileName, outFileName):
    """write elemTable as a python dictionary that can be imported"""

    f = open(outFileName,'w')
    f.write("# File generated from %s\n#\n"%inFileName)
    f.write("babel_elements = {\n")
    for k,v in elemTab.items():
        f.write("  '%s': %s, \n" % (k,str(v)))
    f.write('}\n#END\n');
    f.close()


def read_types_table(filename):
    f = open(filename)
    typestab = {}
    nrow, ncol = map( int, string.split(f.readline()))
    typeFormats = string.split(f.readline())
    for t in typeFormats:
        typestab[t] = []
    for i in range(nrow-1):
        typeNames = string.split(f.readline())
        for j in range(ncol):
            typestab[typeFormats[j]].append(typeNames[j])
    f.close()
    return typestab


def writeTypesTableAsPythonCode(typestab, inFileName, outFileName):
    """write typestab as a python dictionary that can be imported"""

    f = open(outFileName,'w')
    f.write("# File generated from %s\n#\n"%inFileName)
    f.write("babel_types = {\n")
    for k,v in typestab.items():
        f.write("  '%s': %s, \n" % (k,str(v)))
    f.write('}\n#END\n');
    f.close()


if __name__ == '__main__':
    # write tables
    et = read_element_table('element.lis')
    writeElementTableAsPythonCode(et, 'element.lis', 'babelElements.py')

    tt = read_types_table('types.lis')
    writeTypesTableAsPythonCode(tt, 'types.lis', 'babelAtomTypes.py')
