"""
    Utilities for PDB2PQR Suite

    This module provides various utilities for the PDB2PQR suite to be
    imported into other Python scripts.

    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis
"""

__date__ = "23 September 2003"
__author__ = "Todd Dolinsky"
SMALL = 1.0e-7
DIHEDRAL = 57.2958

import string
from math import *

class Matrix:
    """
        Matrix class

        A class for handling matrices
    """
    def __init__(self, lists):
        """
            Create a new matrix object.

            Parameters
                lists:  A list of lists containing the matrix (list)
        """
        self.info = lists
        self.rows = len(lists)
        self.cols = 0
        for row in lists:
            if self.cols == 0: self.cols = len(row)
            if len(row) != self.cols:
                raise ValueError, "Irregularly sized matrix!"

    def __str__(self):
        """
            Print the contents of the matrix

            Returns
                out:  The printed matrix (string)
        """
        out = ""
        for row in self.info:
            for item in row:
                out = "%s %s" % (out, string.rjust(str(item),6))
            out = "%s\n" % out
        return out

    def LU(self, b):
        """
            Solve the matrix Ax = b by LU decomposition:
            Given Ax = b and LU = A,
                Ax = (LU)x = L(Ux) = b
                Solve Ly = b, and then Ux = y. 
            Parameters:
                b = 1 x N matrix, where N is the number of variables
                    (Matrix)
            Returns:
                x = The solved N-element list (list)
        """
        m = self.rows
        n = self.cols
        x = []
        y = []
        
        for i in range(n):
            x.append(0.0)

        # Intialize L to identity, U as a copy of this matrix
        
        ident = []
        for i in range(m):
            list = []
            for j in range(n):
                if i == j: list.append(1.0)
                else: list.append(0.0)
            ident.append(list)
        L = Matrix(ident)
        U = Matrix(self.info)

        # Perform LU decomp

        for i in range(m):
            for j in range(i+1,m):
                if U.info[i][i] == 0.0:
                    raise ValueError, "LU decomposition needs non-zero diags!"
                val = float(U.info[j][i])/U.info[i][i]
                L.info[j][i] = val
                for k in range(n):
                    U.info[j][k] -= U.info[i][k] * val

        # Solve Ly = b, where y = Ux
        
        for i in range(m):
            mult = 1/L.info[i][i]
            sum = b[i]
            for j in range(0,i):
                sum -= L.info[i][j] * y[j]
            y.append((mult*sum))

        # Solve Ux = y
        
        for i in range(m):
            rev = (m-1) - i
            mult = 1/U.info[rev][rev]
            sum = y[rev]
            for j in range(0,i):
                j = (m-1) - j
                sum -= U.info[rev][j] * x[j]
            x[rev] = (mult*sum)
            
        return x

def shortestPath(graph, start, end, path=[]):
    """
        Uses recursion to find the shortest path from one node to
        another in an unweighted graph.  Adapted from
        http://www.python.org/doc/essays/graphs.html .

        Parameters:
            graph: A mapping of the graph to analyze, of the form
                   {0: [1,2], 1:[3,4], ...} . Each key has a list
                   of edges.
            start: The ID of the key to start the analysis from
            end:   The ID of the key to end the analysis
            path:  Optional argument used during the recursive step
                   to keep the current path up to that point

        Returns:
            (variable): Returns a list of the shortest path (list)
                        Returns None if start and end are not
                        connected
    """

    path = path + [start]
    if start == end:
        return path
    if not graph.has_key(start):
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = shortestPath(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest

def analyzeMap(map, value, list=[]):
    """
        Analyze a map of interactions to determine the overall
        connectivity.

        Parameters
            map   : A dictionary of lists which contain the connections
                    to the key (dictionary)
            value : The key value to analyze (variable)

        Returns
            list  : A connectivity list of the map (list)

        Example
            Given map {1: [2], 4: [5], 7: [5,9], 9: [14]} list will return
               For 1:  [1,2]
               For 4,5,7,9,14: [4,5,7,9,14]
               For all other X: [X]
    """
    
    if value in list:
        return []
    else:
        list.append(value)
        
    if value in map:
        for entry in map[value]:
            newlist = analyzeMap(map, entry, list)
            for newitem in newlist:
                if newitem not in list:
                    list.append(newitem)
                    
    for key in map:
        for entry in map[key]:
            if entry == value and key not in list:
                newlist = analyzeMap(map, key, list)
                for newitem in newlist:
                    if newitem not in list:
                        list.append(newitem)
    return list

def getFile(path):
    """
        Obtain a PDB file.  First check the path given on the command
        line - if that file is not available, obtain the file from the
        PDB webserver at http://www.rcsb.org/pdb/ .

        Parameters
            path:  Name of PDB file to obtain (string)

        Returns
            file:  File object containing PDB file (file object)
    """

    import os, urllib

    file = None
    if not os.path.isfile(path):
        URLpath = "http://www.rcsb.org/pdb/cgi/export.cgi/" + path + \
                  ".pdb?format=PDB&pdbId=" + path + "&compression=None"
        file = urllib.urlopen(URLpath)
    else:
        file = open(path)
    return file
        
def distance(coords1, coords2):
    """
        Calculate the distance between two coordinates, as denoted by

            dist = sqrt((x2- x1)^2 + (y2 - y1)^2 + (z2 - z1)^2))

        Parameters
            coords1: Coordinates of form [x,y,z]
            coords2: Coordinates of form [x,y,z]
        Returns
            dist:  Distance between the two coordinates (float)
    """
    dist = 0.0
    list = []

    p = coords2[0] - coords1[0]
    q = coords2[1] - coords1[1]
    r = coords2[2] - coords1[2]
    dist = sqrt(p*p + q*q + r*r)

    return dist

def add(coords1, coords2):
    """
        Add one 3-dimensional point to another
        
        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            list:  List of coordinates equal to coords2 + coords1 (list)
    """
    x = coords1[0] + coords2[0]
    y = coords1[1] + coords2[1]
    z = coords1[2] + coords2[2]
    return [x,y,z]

def subtract(coords1, coords2):
    """
        Subtract one 3-dimensional point from another

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            list:  List of coordinates equal to coords1 - coords2 (list)
    """
    x = coords1[0] - coords2[0]
    y = coords1[1] - coords2[1]
    z = coords1[2] - coords2[2]
    return [x,y,z]

def cross(coords1, coords2):
    """
        Find the cross product of two 3-dimensional points

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            list:  Cross product coords2 and coords1 (list)
    """
    list = []
    x = coords1[1]*coords2[2] -  coords1[2]*coords2[1]
    y = coords1[2]*coords2[0] -  coords1[0]*coords2[2]
    z = coords1[0]*coords2[1] -  coords1[1]*coords2[0]
    list = [x,y,z]
    return list

def dot(coords1, coords2):
    """
        Find the dot product of two 3-dimensional points

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            value:  Dot product coords2 and coords1 (float)
    """
    value = 0.0
    for i in range(3):
        value += coords1[i]*coords2[i]
    return value

def normalize(coords):
    """
        Normalize a set of coordinates
        
        Parameters
            coords: coordinates of form [x,y,z]
        Returns
            list: normalized coordinates (list)
    """
    list = []
    dist = sqrt(pow(coords[0],2) + pow(coords[1],2) + pow(coords[2],2))
    if dist > SMALL:
        a = coords[0]/dist
        b = coords[1]/dist
        c = coords[2]/dist
        list = [a,b,c]
    else:
        list = coords
    return list



def getDihedral(coords1, coords2, coords3, coords4):
    """
        Calculate the angle using the four atoms

        Parameters
            coords1: First of four coordinates of form [x,y,z]
            coords2: Second of four
            coords3: Third of four
            coords4: Fourth of four
        Returns
            value: Size of the angle (float)
    """
    value = 0.0
    
    list43 = subtract(coords4, coords3)
    list32 = subtract(coords3, coords2)
    list12 = subtract(coords1, coords2)

    A = cross(list12, list32)
    Anorm = normalize(A)
    B = cross(list43, list32)
    Bnorm = normalize(B)
    
    scal = dot(Anorm, Bnorm)
    if abs(scal + 1.0) < SMALL:
        value = 180.0
    elif abs(scal - 1.0) < SMALL:
        value = 0.0
    else:
        value = DIHEDRAL * acos(scal)

    chiral = dot(cross(Anorm, Bnorm),list32)
    if chiral < 0:
        value = value * -1.0
    return value
    
def placeOxygen(CA, C, N):
    """
        Place an oxygen according to the planar atoms CA, C, and N using
        a trans-peptide geometry.  Allows for a more accurate method of
        adding oxygen atoms.

        Parameters
            CA:        The coordinates of the CA atom (list)
            C:         The coordinates of the C atom (list)
            N:         The coordinates of the peptide bonded N atom from the
                       next residue (list)
        Returns
            location:  The location of the residue (list)
    """

    # Step 1: Find the vector normal to the C-CA-N plane in order to get
    #         the equation for any point in the plane

    vec1 = subtract(CA,C)
    vec2 = subtract(N,C)
    planeeq = cross(vec1, vec2)
    sum = planeeq[0]*C[0] + planeeq[1]*C[1] + planeeq[2]*C[2]

    # Step 2: Get two more equations using the known C-O distance (1.24 A)
    #         and CA-C-O and N-C-O bond angles (120.5 and 123.5 degrees,
    #         respectively) using the identity
    #
    #         A . B = |A||B| cos(angle)
    
    num1 = math.sqrt(pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2)) 
    num2 = math.sqrt(pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2))
    
    # For vector 1
    
    val1 = 0
    angle = 120.5*math.pi/180.0
    val1 = num1*1.24*math.cos(angle)
    for j in range(3):
        val1 += C[j]*vec1[j]
        
    # For vector 2

    val2 = 0
    angle = 123.5*math.pi/180.0
    val2 = num2*1.24*math.cos(angle)
    for j in range(3):
        val2 += C[j]*vec2[j]
        
    # Step 3: We now use Gaussian Elimination to solve the following matrix
    #
    #         [ planeq[0] planeeq[1] planeeq[2] ]  =  [sum]
    #         [  vec1[0]    vec1[1]  vec1[2]    ]  =  [val1]
    #         [  vec2[0]    vec2[1]  vec2[2]    ]  =  [val2]
    

    fac1 = -1 * planeeq[0]/vec1[0]
    new1 = [0, fac1*vec1[1]+planeeq[1], fac1*vec1[2]+planeeq[2]]
    val1 = fac1*val1+sum
    
    fac2 = -1 * planeeq[0]/vec2[0]
    new2 = [0, fac2*vec2[1]+planeeq[1], fac2*vec2[2]+planeeq[2]]
    val2 = fac2*val2+sum
    
    fac3 = -1 * new1[1]/new2[1]
    newest = [0,0,fac3*new2[2]+new1[2]]
    val2 = fac3*val2+val1

    # Step 4: Backfill in to find the results
    
    z = val2/newest[2]
    y = (val1 - z*new1[2])/new1[1]
    x = (sum - z*planeeq[2] - y*planeeq[1])/planeeq[0]

    location = [x,y,z]
    return location
