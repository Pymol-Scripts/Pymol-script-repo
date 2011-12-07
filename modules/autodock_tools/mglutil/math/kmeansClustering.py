import random
import numpy
from math import sqrt
# -- The Point class represents points in n-dimensional space

class Point:
    # Instance variables
    # self.coords is a list of coordinates for this Point
    # self.n is the number of dimensions this Point lives in (ie, its space)
    # self.reference is an object bound to this Point
    # Initialize new Points

    def __init__(self, coords, reference=None):
        self.coords = coords
        self.n = len(coords)
        self.reference = reference
    # Return a string representation of this Point

    def __repr__(self):
        return str(self.coords)
# -- The Cluster class represents clusters of points in n-dimensional space


class Cluster:
    # Instance variables
    # self.points is a list of Points associated with this Cluster
    # self.n is the number of dimensions this Cluster's Points live in
    # self.centroid is the sample mean Point of this Cluster

    def __init__(self, points):
        # We forbid empty Clusters (they don't make mathematical sense!)
        if len(points) == 0: raise Exception("ILLEGAL: EMPTY CLUSTER")
        self.points = points
        self.n = points[0].n
        # We also forbid Clusters containing Points in different spaces
        # Ie, no Clusters with 2D Points and 3D Points
        for p in points:
            if p.n != self.n: raise Exception("ILLEGAL: MULTISPACE CLUSTER")
        # Figure out what the centroid of this Cluster should be
        self.centroid = self.calculateCentroid()
    # Return a string representation of this Cluster

    def __repr__(self):
        return str(self.points)
    # Update function for the K-means algorithm
    # Assigns a new list of Points to this Cluster, returns centroid difference

    def update(self, points):
        old_centroid = self.centroid
        self.points = points
        self.centroid = self.calculateCentroid()
        x1,y1,z1 = old_centroid.coords
        x2,y2,z2 = self.centroid.coords
        return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) )

    # Calculates the centroid Point - the centroid is the sample mean Point
    # (in plain English, the average of all the Points in the Cluster)
    def calculateCentroid(self):
        centroid_coords = []
        # For each coordinate:
        for i in range(self.n):
            # Take the average across all Points
            centroid_coords.append(0.0)
            for p in self.points:
                centroid_coords[i] = centroid_coords[i]+p.coords[i]
            centroid_coords[i] = centroid_coords[i]/len(self.points)
        # Return a Point object using the average coordinates
        return Point(centroid_coords)

    def radiusOfGyration(self):
        ptCoords = [x.coords for x in self.points]
        delta = numpy.array(ptCoords)-self.centroid.coords
        rg = sqrt( sum( numpy.sum( delta*delta, 1))/float(len(ptCoords)) )
        return rg

    def encapsualtingRadius(self):
        ptCoords = [x.coords for x in self.points]
        delta = numpy.array(ptCoords)-self.centroid.coords
        rM = sqrt( max( numpy.sum( delta*delta, 1)) )
        return rM


# -- Return Clusters of Points formed by K-means clustering
def kmeans(points, k, cutoff, initial=None):
    # Randomly sample k Points from the points list, build Clusters around them
    if initial is None:
        # Randomly sample k Points from the points list, build Clusters around them
        initial = random.sample(points, k)
    else:
        assert len(initial)==k

    clusters = []
    for p in initial: clusters.append(Cluster([p]))
    # Enter the program loop
    while True:
        # Make a list for each Cluster
        lists = []
        for c in clusters: lists.append([])
        # For each Point:
        for p in points:
            # Figure out which Cluster's centroid is the nearest
            x1,y1,z1 = p.coords
            x2,y2,z2 = clusters[0].centroid.coords
            smallest_distance = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                                      (z1-z2)*(z1-z2) )
            index = 0
            for i in range(len(clusters[1:])):
                x2,y2,z2 = clusters[i+1].centroid.coords
                distance = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                                 (z1-z2)*(z1-z2) )
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i+1
            # Add this Point to that Cluster's corresponding list
            lists[index].append(p)
        # Update each Cluster with the corresponding list
        # Record the biggest centroid shift for any Cluster
        biggest_shift = 0.0
        for i in range(len(clusters)):
            if len(lists[i]):
                shift = clusters[i].update(lists[i])
                biggest_shift = max(biggest_shift, shift)
        # If the biggest centroid shift is less than the cutoff, stop
        if biggest_shift < cutoff: break
    # Return the list of Clusters
    return clusters
