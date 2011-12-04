## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
# $Id: AutoDockBondClassifier.py,v 1.9 2008/08/13 19:08:24 rhuey Exp $
#
#
#
#


from MolKit.bondClassifier import BondClassifier
from MolKit.bondSelector import AmideBondSelector, PeptideBackBoneBondSelector
from MolKit.bondSelector import HydrogenRotatorBondSelector, RotatableBondSelector
from MolKit.bondSelector import LeafBondSelector, CycleBondSelector
from MolKit.bondSelector import GuanidiniumBondSelector, BondOrderBondSelector 
from MolKit.bondSelector import AromaticCycleBondSelector2
from MolKit.molecule import BondSet
       
def dist(c1, c2):       
    import numpy.oldnumeric as Numeric, math
    d = Numeric.array(c2) - Numeric.array(c1) 
    ans = math.sqrt(Numeric.sum(d*d))
    return round(ans, 3)


class AutoDockBondClassifier(BondClassifier):

    def __init__(self, tolerance=0.01, detectAll=True):
        self.detect_all_cycles = detectAll
        self.d = {
            'amide':AmideBondSelector(),
            'ppbb':PeptideBackBoneBondSelector(),
            'leaf':LeafBondSelector(),
            'cycle':CycleBondSelector(),
            'rotatable':RotatableBondSelector(),
            'bondOrder2':BondOrderBondSelector(2),
            'hydrogenRotators':HydrogenRotatorBondSelector(),
            'guanidinium':GuanidiniumBondSelector(),
            'aromatic': AromaticCycleBondSelector2()

            }
        BondClassifier.__init__(self, self.d)
        #used to detect colinear atoms
        #if dist1+dist2<dist13+0.1
        self.tolerance = 0.01
#NOTE hydrogenRotators value is a group of ATOMS not BONDS


    def classify(self, bonds=None):
        """ 
        select using each bondselector (the values of the dict); store result
        in resultDict and return the result dict when finished...
        """
        #make sure that classify is called with some bonds
        assert isinstance(bonds, BondSet)
        resultDict = {}
        rotatables = self.d['rotatable'].select(bonds)

        for k, v in self.dict.items():
            if k=='bondOrder2':
                resultDict[k] = v.select(bonds,2)
            elif k=='cycle':
                resultDict[k] = v.select(bonds,self.detect_all_cycles)
            else:
                resultDict[k] = v.select(bonds)
        bnds = resultDict['leaf']
        for b in bnds:
            at1 = b.atom1
            at2 = b.atom2
            if len(at1.bonds)==1:
                at1.leaf = 1
                if len(at2.bonds)==2:
                    dist1 = dist(at1.coords, at2.coords)
                    #find neighbor of at2 and check bond distances
                    b2 = at2.bonds[0]
                    if b2==b:
                        b2 = at2.bonds[1]
                    neighbor = b2.atom1
                    if neighbor==at2:
                        neighbor = b2.atom2
                    dist2 = dist(neighbor.coords, at2.coords)
                    dist13 = dist(at1.coords, neighbor.coords)
                    if dist13+self.tolerance>=dist1+dist2:
                        #print "adding ", b2, 'to leaf bonds'
                        resultDict['leaf'].append(b2)
                        if b2 in resultDict['rotatable']:
                            resultDict['rotatable'].remove(b2)
            elif len(at2.bonds)==1:
                at2.leaf = 1
                if len(at1.bonds)==2:
                    dist1 = dist(at1.coords, at2.coords)
                    #find neighbor of at1 and check bond distances
                    b2 = at1.bonds[0]
                    if b2==b:
                        b2 = at1.bonds[1]
                    neighbor = b2.atom1
                    if neighbor==at1:
                        neighbor = b2.atom2
                    dist2 = dist(neighbor.coords, at1.coords)
                    dist13 = dist(at2.coords, neighbor.coords)
                    if dist13+self.tolerance>=dist1+dist2:
                        resultDict['leaf'].append(b2)
                        #print "2: adding ", b2, 'to leaf bonds'
                        b2.possibleTors=0
                        if b2 in resultDict['rotatable']:
                            resultDict['rotatable'].remove(b2)
                    
        #8/16:NOTE only select from rotatables for amide+ guanidinium 
        for k in ['amide', 'guanidinium']:
            resultDict[k] = self.d[k].select(rotatables)
        return resultDict


