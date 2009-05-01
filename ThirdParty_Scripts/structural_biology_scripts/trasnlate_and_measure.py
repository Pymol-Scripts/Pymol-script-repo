def translateAndMeasure(selection, other, translationVector, cutoff):
    cmd.translate(translationVector, selection)
    return checkDistances(selection, other, cutoff)
 
def checkDistances(moleculeA, moleculeB, cutoff):
    ids_A = getIds(moleculeA)
    ids_B = getIds(moleculeB)
    for idA in ids_A:
        for idB in idsB:
            d = distance(moleculeA, idA, moleculeB, idB)
            if d > cutoff: return "overlap"
    return "no overlap"
 
def distance(a, idA, b, idB):
    atomA = "%s and id %s" % (a, idA)
    atomB = "%s and id %s" % (b, idB)
    return cmd.get_distance(atomA, atomB)
 
def getIds(selection):
    my_dict = { 'my_list' : [] }
    cmd.iterate(selection, "my_list.append(ID)", space=my_dict)
    return my_dict['my_list']

