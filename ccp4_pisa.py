'''
See more here: http://www.pymolwiki.org/index.php/ccp4_pisa

        parse a PISA contact file and create atoms selections for all interfaces

        For each interface, two selections are created containing the atoms of
        the interface on each chain. The selection names follow the convention
        bondtype_#Id_Chain1Chain2. If Chain1 equals Chain2 then the two selections are
        numbered.

        bondtype corresponds to h-bonds, salt-bridges, ss-bonds and cov-bonds and
        are marked with the prefixes hb, sb, ss and cov.

        For example, all h-bonds in interface 3 between chain A and D create the selections
        hb_3_AD and hb_3_DA.

        Salt-bridges in interface 4 between chain A and a symmetry copy of A creates the selections
        sb_4_AA1 and sb_4_AA2.

    PARAMS
        filename
            filename of the PISA contacts file

    RETURNS
        a set of selections in PYMOL.

    AUTHOR
        Gerhard Reitmayr and Dalia Daujotyte, 2011.
'''
from __future__ import print_function
from pymol import cmd
from xml.etree import ElementTree


def parseElement(element):
    """ creates a dict for the sub elements of the element"""
    result = {}
    for l in range(len(element)):
        if element[l].text != None:
            result[element[l].tag.strip()] = element[l].text.strip()
    return result


def parseBond(elementDir):
    """ puts bond information into tuples"""
    return ((elementDir['chain-1'], elementDir['seqnum-1'], elementDir['atname-1']), (elementDir['chain-2'], elementDir['seqnum-2'], elementDir['atname-2']))


def parseInterface(interface, bondname):
    """ parses a single interface into the interface id, the two chain names connected
        and two lists of atoms for each chain"""
    bonds = interface.findall(bondname)
    id = interface.find('id').text.strip()
    if len(bonds) == 0:
        return None
    left = []
    right = []
    for b in bonds:
        l, r = parseBond(parseElement(b))
        left.append(l)
        right.append(r)
    left_chain = left[0][0]
    right_chain = right[0][0]
    if left_chain > right_chain:
        return id, (right_chain, left_chain), right, left
    return id, (left_chain, right_chain), left, right


def createSelectionList(atomlist):
    """creates a PYMOL selection string for a list of atoms"""
    atomnames = [chain + '/' + res + '/' + atom for chain, res, atom in atomlist]
    return " or ".join(atomnames)


def createInterfaceSelection(interface, prefix):
    """creates two selections for an interfaces"""
    id, e, l, r = interface
    leftname = prefix + '_' + str(id) + '_' + e[0] + e[1]
    rightname = prefix + '_' + str(id) + '_' + e[1] + e[0]
    if e[0] == e[1]:
        leftname = leftname + '1'
        rightname = rightname + '2'
    leftlist = createSelectionList(l)
    rightlist = createSelectionList(r)
    try:
        cmd.select(leftname, leftlist)
        cmd.select(rightname, rightlist)
    except:
        print(leftname, '\t', leftlist)
        print(rightname, '\t', rightlist)
    return leftname, rightname


def ccp4_pisa(filename):
    bond_types = [
        ('h-bonds', 'h-bonds/bond', 'hb'),
        ('salt-bridges', 'salt-bridges/bond', 'sb'),
        ('ss-bonds', 'ss-bonds/bond', 'ss'),
        ('cov-bonds', 'cov-bonds/bond', 'cv')
    ]

    tree = ElementTree.parse(open(filename))
    interfaces = tree.findall('//interface')

    result = []
    for name, path, prefix in bond_types:
        allcontacts = [edge for edge in [parseInterface(i, path) for i in interfaces] if edge != None]

        allselections = []
        for c in allcontacts:
            allselections.extend(createInterfaceSelection(c, prefix))

        result.append(len(allselections) / 2)
        if len(allselections) > 0:
            try:
                cmd.select(name, " or ".join(allselections))
            except:
                print(name, '\t', " or ".join(allselections))

    print('selectPISAContacts found interfaces with', end=' ')
    for number, type in zip(result, bond_types):
        print(number, type[0], ",", end=' ')

try:
    cmd.extend("ccp4_pisa", ccp4_pisa)
except:
    # for debugging
    ccp4_pisa('../pisa/interfaces_2c7r.pisa')
