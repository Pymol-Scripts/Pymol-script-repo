#########################################################################
#
# Date: Jan 2004  Author: Daniel Stoffler
#
# stoffler@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler and TSRI
#
#########################################################################

import os, sys, string
from DejaVu.IndexedPolygons import IndexedPolygons
from mglutil.util.packageFilePath import which

def findQhullModule(modname):
    """return (platform-dependent) path to qhull module"""

    if modname not in ['qconvex', 'qdelaunay', 'qhalf', 'qhull',
                       'qvoronoi', 'rbox']:
        print 'QHULL ERROR! Illegal module name %s.'%modname
        return None
    try:
        mod = __import__('binaries')
        pth = mod.__path__[0]
    except ImportError:
        pth = ''

    pth = os.path.join( pth, modname)

    if which(pth):
        #print pth
        return pth
    else:
        #print 'QHULL ERROR! Module not found in %s'%pth
        return None
    

class QConvex:
    """Compute convex hull based on a list of 3-D coordinates. Return a
DejaVu IndexedPolygon geometry.

Usage: self.computeConvex(coords, tpath)
  coords: a list of [x,y,z] values
  tpath = optional path for temporary files
"""


    def __init__(self, name='qconvex'):
        self.geom = IndexedPolygons(name, inheritMaterial=0)
        self.tmpPath = './' # user can specify path for tmp files
                            # using setTmpPath() method 
        

    #------------------ HELPER FUNCTIONS ---------------------------------#
 
    def writeQConvex(self, filename, coords):
        """QConvex: http://www.qhull.org/html/qconvex.htm
Input a filename and [x,y,z]-coordinates, save this in QConvex format."""

        data = []
        data.append('3 RBOX c\n')
        data.append( str(len(coords))+'\n')
        for (x,y,z) in coords:
            data.append( '%f %f %f\n'%(x,y,z))

        try:
            f = open( os.path.join(self.tmpPath, filename), 'w')
            f.writelines(data)
            f.close()
        except:
            print 'QCONVEX ERROR! Cannot write into %s'%self.tmpPath


    def readQConvex(self, filename):
        """QConvex: http://www.qhull.org/html/qconvex.htm

Read a QConvex output file, output a DejaVu IndexedPolygon.

 Data Format:
    [...]print vertices and facets of the convex hull in OFF format. The first
line is the dimension. The second line is the number of vertices, facets, and
ridges. The vertex coordinates are next, followed by the facets. Each facet
starts with the number of vertices. The cube example has four vertices per
facet."""
        
        try:
            f = open( os.path.join(self.tmpPath, filename), 'r')
            data = f.readlines()
            f.close()
        except:
            print 'QCONVEX ERROR! Temp. file not found in %s'%self.tmpPath
            return
        
        # get more header info
        header = string.split(data[1])
        lenVerts = int(header[0])
        lenFaces = int(header[1])

        vertices = []
        faces = []

        for d in data[2:lenVerts+2]: # offset of 2 because of file header
            spl = string.split(d)
            vertices.append( [float(spl[0]), float(spl[1]), float(spl[2])] )

        for d in data[lenVerts+2:]: # offset of 2 because of file header
            spl = map( int, string.split(d) )

            for i in range(3, len(spl)):
                faces.append( [spl[1], spl[i], spl[i-1]] )

        self.geom.Set(vertices=vertices, faces=faces)
        

    def setTmpPath(self, path):
        """set the path for the two temporary files. Note: if the specified
        path does not exist, we try to write into the startup directory"""
        
        if path is None:
            path = './'
        if path == '':
            path = os.path.abspath("./")
            
        if not os.path.exists(path):
            print 'QCONVEX ERROR! Path %s does not exist!'%path
            # use path where we started Python process
            self.tmpPath = os.path.join( os.path.abspath('./'), '')
            print 'Trying to save temp. file in: %s'%self.tmpPath

        else:
            self.tmpPath = os.path.join( os.path.abspath(path), '')


    def getGeom(self):
        """returns the DejaVu IndexedPolygon geometry"""
        return self.geom

    #------------------ END HELPER FUNCTIONS ------------------------------#


    #------------------ USE THIS METHOD: ----------------------------------#
    def computeConvex(self, coords, tpath=None):

        ### 1) clean directory, save file in qconvex format:
        if tpath is not None:
            self.setTmpPath(tpath)
        try:
            os.remove( os.path.join(self.tmpPath, 'tmp_qconvex_input') )
        except:
            pass
        try:
            os.remove( os.path.join(self.tmpPath, 'tmp_qconvex_output') )
        except:
            pass

        self.writeQConvex(
            os.path.join(self.tmpPath, 'tmp_qconvex_input'), coords)

        ### 2) run qconvex, create new output file

        # FIXME: Please note, this is a hack: we build the path where
        # qconvex resides, depending on the operating system. In the future,
        # qconvex might be accessed through a SOAP service!

        # build a path string for the file qconvex
        pth = findQhullModule('qconvex')
        if pth is None:
            return

        # build string to be executed
        execstring = pth + ' o < ' + self.tmpPath + 'tmp_qconvex_input > '+\
        self.tmpPath + 'tmp_qconvex_output'
        os.system(execstring)

        ### 3) load output, return IndexedPolygon
        self.readQConvex(
            os.path.join(self.tmpPath, 'tmp_qconvex_output') )

        ### 4) clean up temporary files
        try:
            os.remove( os.path.join(self.tmpPath, 'tmp_qconvex_input') )
        except:
            print 'Cannot delete temporary input file'
        try:
            os.remove( os.path.join(self.tmpPath, 'tmp_qconvex_output') )
        except:
            print 'Cannot delete temporary output file'



        
