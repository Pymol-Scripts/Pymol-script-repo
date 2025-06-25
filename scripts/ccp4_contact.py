'''
See more here: http://www.pymolwiki.org/index.php/ccp4_contact

    ccp4_contact -- parses CCP4/CONTACT log file and selects residues and atoms.
    http://www.ccp4.ac.uk/html/contact.html
 
    PARAMS
        contactsfile
            filename of the CCP4/CONTACT log file
 
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
        Gerhard Reitmayr and Dalia Daujotyte, 2011.
'''
from pathlib import Path
from pymol import cmd
import re


def parseCONTACTContacts(f):
    # Lys    24A  ca  Asp   263D  CG   ...  4.94    [   -1B   ]   3: -X,  Y+1/2,  -Z+1/2
    conParser = re.compile("(\S*)\s*(\d+)([A-Z])\s*(\w+)")
    s1 = []
    s2 = []
    for line in f:
        if line.startswith('#'):
            continue
        matches = conParser.findall(line)
        if len(matches) == 2:
            s1.append(matches[0])
            s2.append(matches[1])
        elif len(matches) == 1:
            s2.append(matches[0])

    return (s1, s2)


def ccp4_contact(contactsfile, selName1="source", selName2="target"):
    # read and parse contacts file into two lists of contact atoms and contact pair list
    s1, s2 = parseCONTACTContacts(open(contactsfile))

    # create a selection for the first contact list

    # create the PYMOL selection macros for the residues
    resNames = [chain + "/" + residue + "/" for (type, residue, chain, atom) in s1]
    # put them in a set to remove duplicates and then join with 'or'
    resSel = " or ".join(frozenset(resNames))
    # finally select them under the new name
    cmd.select(selName1 + "_res", resSel)

    atomNames = [chain + "/" + residue + "/" + atom for (type, residue, chain, atom) in s1]
    atomSel = " or ".join(frozenset(atomNames))
    cmd.select(selName1 + "_atom", atomSel)

    # create a selection for the second contact list

    resNames = [chain + "/" + residue + "/" for (type, residue, chain, atom) in s2]
    resSel = " or ".join(frozenset(resNames))
    cmd.select(selName2 + "_res", resSel)

    atomNames = [chain + "/" + residue + "/" + atom for (type, residue, chain, atom) in s2]
    atomSel = " or ".join(frozenset(atomNames))
    cmd.select(selName2 + "_atom", atomSel)

cmd.extend("ccp4_contact", ccp4_contact)


def test_ccp4_contact():
    datadir = Path(__file__).parents[1] / "files_for_examples"
    cmd.reinitialize()
    cmd.fab("PCQAFSISGKQKGFEDSRGTL", chain="A", resi=80)
    ccp4_contact(f"{datadir}/2c7r.contact", selName1="sel1", selName2="sel2")
    assert cmd.count_atoms("sel1_res") == 128
    assert cmd.count_atoms("sel1_atom") == 20
    assert cmd.count_atoms("sel2_res") == 128
    assert cmd.count_atoms("sel2_atom") == 20