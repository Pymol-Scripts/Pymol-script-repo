import cmd
from chempy import models, cpv
 
"""
Create an atom at a distance 'distance' along the bond between atomA and atomB
"""
def createAtomAlongBond(modelName, distance, resiA, atomNameA, resiB, atomNameB, atomNameC):
    model = cmd.get_model(modelName)
    p1 = getAtomCoords(model, str(resiA), atomNameA)
    p2 = getAtomCoords(model, str(resiB), atomNameB)
    if p1 is None:
        print "atom not found!", modelName, resiA, atomNameA
    elif p2 is None:
        print "atom not found!", modelName, resiB, atomNameB
    else:
        p3 = calculateNewPoint(p1, p2, distance)
 
        # the details of the new atom
        atomDetails = {}
        atomDetails['residueName'] = "HOH"
        atomDetails['residueNumber'] = "1"
        atomDetails['symbol'] = "O"
        atomDetails['name'] = atomNameC
        atomDetails['coords'] = p3
 
        # make an atom with index n+1 and chain "X"
        newAtom = makeAtom(model.nAtom + 1, atomDetails, "X")
        model.add_atom(newAtom)
        model.update_index()
        cmd.load_model(model, "newpeptide")
 
def getAtomCoords(model, resi, atomName):
    for a in model.atom:
        if a.resi == resi and a.name == atomName:
            return a.coord
    return None
 
def calculateNewPoint(p1, p2, distance):
    v1 = cpv.normalize(cpv.sub(p1, p2))
    return cpv.add(p1, cpv.scale(v1, distance))
 
def makeAtom(index, atomDetails, chain):
    atom = chempy.Atom()
    atom.index = index
    atom.name = atomDetails['name']
    atom.symbol = atomDetails['symbol']
    atom.resn = atomDetails['residueName']
    atom.chain = chain
    atom.resi = atomDetails['residueNumber']
    atom.resi_number = int(atomDetails['residueNumber'])
    atom.coord = atomDetails['coords']
    atom.hetatm = False
    return atom
 
cmd.extend("createAtomAlongBond", createAtomAlongBond)
