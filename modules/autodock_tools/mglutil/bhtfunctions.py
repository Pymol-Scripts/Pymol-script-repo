## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

import numpy.oldnumeric as Numeric
from bhtree import bhtreelib
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/bhtfunctions.py,v 1.9 2007/07/24 17:30:40 vareille Exp $
#
# $Id: bhtfunctions.py,v 1.9 2007/07/24 17:30:40 vareille Exp $
#

# ClosePointsDist2: result and dist are empty arrays large enough to contain
# all the points expected to be found. To be safe, they should be the same
# size as the list of coordinates. This function then puts in the results
# array the indices of the close points in the list (as supplied in the array
# ids): the dist array contains the corresponding distances.

def findNearestAtoms(mol,vertices, **kw):
    """None <- color(mol,vertices2,**kw)
    mol:           reference molecule
    vertices:      list of lists(coordinates): the first three items in each list
                   must be coordinates x,y,z of a point.
                   
    atomIndices is the index of the nearest atom to the vertex, such that
    mol.allAtoms[atomIndices[x]] is the nearest atom to vertices[x]
    vertexIndices is the list of nearest vertices to an atom, such that
    vertexIndices[x] = [vertex1,vertex2,...] are the vertices associated with
    mol.allAtoms[x]
    """
    
    coords = mol.allAtoms.coords
    if not hasattr(mol,'bhtree'):
        print "Building bhtree for ",mol
        ids = Numeric.arange(len(coords)).astype('i')
        bhtree = bhtreelib.TBHTree(coords,ids,10,10,9999.0)
        mol.bhtree = bhtree
    
    vertexIndices={}
    atomIndices={}
    for x in range(len(coords)):
        vertexIndices[x+1]=[]

    cutoff=5.
    for x in range(len(vertices)):
        xyz = vertices[x]
        result = Numeric.zeros( (len(vertices),) ).astype('i')
        dist = Numeric.zeros( (len(vertices),) ).astype('f')
        nb2 = mol.bhtree.ClosePointsDist2(tuple(xyz[:3]), cutoff, result, dist )
        while nb2==0:
            cutoff = cutoff+5.
            nb2 = mol.bhtree.ClosePointsDist2(tuple(xyz[:3]), cutoff, result, dist )
        result = result[:nb2]
        dist = dist[:nb2]
        idx = dist.tolist().index(min(dist))
        atnum = result[idx]+1
        atomIndices[x]=atnum
        vertexIndices[atnum].append(x)
        
    return atomIndices,vertexIndices








