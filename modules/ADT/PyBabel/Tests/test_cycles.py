import sys, string
from mglutil.regression import testplus

class Atom:

    def __init__(self, coords):
        self.coords = coords
        self.bonds = []

class Bond:

    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom1.bonds.append(self)
        self.atom2.bonds.append(self)


def getGraphFromFile(filename):
    f = open(filename)
    nba = int(f.readline())
    atoms = []
    for i in range(nba):
        data = string.split(f.readline())
        x,y,z = map( float, data)
        atoms.append(Atom((x,y,z)))
    bonds = []
    nbb = int(f.readline())
    for i in range(nbb):
        data = string.split(f.readline())
        at1,at2 = map( int, data)
        bonds.append( Bond( atoms[at1], atoms[at2] ) )
    return atoms, bonds


def test_findCycles4():
    # test that atoms in cycles are ordered on taxol 
    atoms, bonds = getGraphFromFile('Data/txp.graph')
    from PyBabel.cycle import RingFinder
    rings = RingFinder()
    rings.findRings2(atoms, bonds)

    def hasBond(a1, a2, cycleHasBond):
        for b in cyclebonds:
            if b.atom1==a1 and b.atom2==a2:
                return 1
            elif b.atom1==a2 and b.atom2==a1:
                return 1
        return 0
    
    def isOrdered(cycleatoms, cyclebonds):

        def hasBond(a1, a2, cyclebonds):
            for b in cyclebonds:
                if b.atom1==a1 and b.atom2==a2:
                    return 1
                elif b.atom1==a2 and b.atom2==a1:
                    return 1
            return 0

        start = cycleatoms[0]
        for end in cycleatoms[1:]:
            if not hasBond(start, end, cyclebonds):
                return 0
            start = end
        return 1
    
    for r in rings.rings:
        assert isOrdered(r['atoms'], r['bonds'])


def test_findCycles1():
    # test cycle detection on taxol 
    atoms, bonds = getGraphFromFile('Data/txp.graph')
    from PyBabel.cycle import RingFinder
    rings = RingFinder()
    rings.findRings2(atoms, bonds)
    assert rings.ringCount==6
    
def test_findCycles2():
    # test cycle detection on modified taxol to have 5 membered ring
    atoms, bonds = getGraphFromFile('Data/fuzedRings2.graph')
    from PyBabel.cycle import RingFinder
    rings = RingFinder()
    rings.findRings2(atoms, bonds)
    assert rings.ringCount==4

    
def test_findCycles3():
    # test cycle detection on taxol merged rings only
    atoms, bonds = getGraphFromFile('Data/fuzedRings1.graph')
    from PyBabel.cycle import RingFinder
    rings = RingFinder()
    rings.findRings2(atoms, bonds)
    assert rings.ringCount==4


harness = testplus.TestHarness( "PyBabel_cycles",
                                funs = testplus.testcollect( globals()),
                                )

if __name__ == '__main__':
    testplus.chdir()
    print harness
    sys.exit( len( harness))
