'''
See more here: http://www.pymolwiki.org/index.php/ccp4_ncont

    ccp4_ncont -- parses CCP4/NCONT log file and selects residues and atoms.
    http://www.ccp4.ac.uk/html/ncont.html

    PARAMS
        contactsfile
            filename of the CCP4/NCONT contacts log file

        selName1
            the name prefix for the _res and _atom selections returned for the
            source set of chain

        selName2
            the name prefix for the _res and _atom selections returned for the
            target set of chain

    RETURNS
        4 selections of interface residues and atoms are created and named
        depending on what you passed into selName1 and selName2

    AUTHOR
        Gerhard Reitmayr and Dalia Daujotyte, 2009.
'''
from __future__ import print_function
from pymol import cmd
import re


def parseNCONTContacts(f):
    # /1/B/ 282(PHE). / CE1[ C]:  /1/E/ 706(GLN). / O  [ O]:   3.32
    # * in the second group is needed when chain code is blank
    conParser = re.compile("\s*/(\d+)/([a-zA-Z0-9]*)/\s*(\d+).*?/\s*([a-zA-Z0-9]*).*?:")
    mode = 0
    s1 = []
    s2 = []
    pairs = []
    for line in f:
        if mode == 0:
            if line.strip().startswith("SOURCE ATOMS"):
                mode = 1
        elif mode == 1:
            mode = 2
        elif mode == 2:
            matches = conParser.findall(line)
            if len(matches) == 0:
                return (s1, s2, pairs)
            if len(matches) == 2:
                s1.append(matches[0])
                s2.append(matches[1])
            elif len(matches) == 1:
                s2.append(matches[0])
            pairs.append((len(s1) - 1, len(s2) - 1))
        else:
            print("Unknown mode", mode)


def ccp4_ncont(contactsfile, selName1="source", selName2="target"):
    # read and parse contacts file into two lists of contact atoms and contact pair list
    s1, s2, pairs = parseNCONTContacts(open(contactsfile))
    # create a selection for the first contact list

    # create the PYMOL selection macros for the residues
    resNames = [chain + "/" + residue + "/" for (type, chain, residue, atom) in s1]
    # put them in a set to remove duplicates and then join with 'or'
    resSel = " or ".join(frozenset(resNames))
    # finally select them under the new name
    cmd.select(selName1 + "_res", resSel)

    atomNames = [chain + "/" + residue + "/" + atom for (type, chain, residue, atom) in s1]
    atomSel = " or ".join(frozenset(atomNames))
    cmd.select(selName1 + "_atom", atomSel)

    # create a selection for the second contact list

    resNames = [chain + "/" + residue + "/" for (type, chain, residue, atom) in s2]
    resSel = " or ".join(frozenset(resNames))
    cmd.select(selName2 + "_res", resSel)

    atomNames = [chain + "/" + residue + "/" + atom for (type, chain, residue, atom) in s2]
    atomSel = " or ".join(frozenset(atomNames))
    cmd.select(selName2 + "_atom", atomSel)

cmd.extend("ccp4_ncont", ccp4_ncont)
