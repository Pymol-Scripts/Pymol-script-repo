import math
from pymol import querying
from pymol.cgo import *
from pymol import cmd
 
#NOTE!! : These packages (numarray, Scientific) must be present in pymol's 
#'$PYMOLDIR/py23/lib/python2.3/site-packages/' directory!!
# OR..if you are using Mac PyMOL 0.99beta19 or later then the system installs will be used
from numarray import *
from numarray.ma import average
from numarray import linear_algebra as la
 
from Scientific.Geometry import Vector, Tensor, Transformation
 
def boundingBox(selectionName, colourRGB=[1,1,1]):
        """
        The main function to call : 
 
                run "box.py"
                boundingBox("peptide")
 
        Should make a box around "peptide" (assuming it exists!). For a different colour use:
 
                boundingBox("peptide", colourRGB=[1, 0, 1])
 
        Or whatever. The box should be a cgo called 'peptide-box' (for this example).
        """
        model = querying.get_model(selectionName)
        coords = model.get_coord_list()
 
        #find the least square plane through the coordinates
        eigenvalues, eigenvectors, centroid = findPlaneWithEigenVectors(coords)
        normal = eigenvectors[eigenvalues.argmin()]
        eigenvectors.remove(normal)
 
        #orient the axes correctly
        x, y, normal = orientAxes(eigenvectors[0], eigenvectors[1], normal)
 
        #determine the dimensions and the structure's orientation wrt the origin
        minDimensions, rotation = findMinDimensionsAndRotation(coords, centroid, x, y, normal)
 
        #'create' the box(IE: make the corners) and 'draw' it (IE: make a cgo)
        box = makeBox(minDimensions, rotation, centroid)
        drawBox(box, selectionName, colourRGB)
 
def findPlaneWithEigenVectors(coords):
        centroid = average(coords)
        coords -= centroid
        B = matrixmultiply(transpose(coords), coords)
        eigenvalues, eigenvectors = la.eigenvectors(B)
        #return eigenvalues, [Vector(e) for e in eigenvectors], Vector(centroid)
        return eigenvalues, [Vector([i for i in e]) for e in eigenvectors], Vector(centroid) #not sure why I had to make this change!
 
def orientAxes(x, y, z):
        XcrossY = x.cross(y)
        #ZXY = around(math.degrees(z.angle(XcrossY)))
        ZXY = int(around(math.degrees(z.angle(XcrossY)))) #again, a bit of a hack!
        if (ZXY == 180): x, y = y, x
        return x, y, z
 
def makeBox(dimensions, rotation, centroid):
        x, y, z = dimensions
        v = [[]] * 8
 
        #make a cuboid with the lower corner on the origin
        v[0] = [0, 0, 0]                # [0, 0, 0]
        v[1] = [x, 0, 0]                # [1, 0, 0]
        v[2] = [x, y, 0]                # [1, 1, 0]
        v[3] = [x, 0, z]                # [1, 0, 1]
        v[4] = [0, y, 0]                # [0, 1, 0]
        v[5] = [0, 0, z]                # [0, 0, 1]
        v[6] = [0, y, z]                # [0, 1, 1]
        v[7] = [x, y, z]                # [1, 1, 1]
 
        #move to the origin, THEN move to the centroid of the points, then rotate
        translationToOrigin = Transformation.Translation(-Vector(x/2, y/2, z/2))
        translationToCentroid = Transformation.Translation(centroid)
        transform = translationToCentroid * rotation * translationToOrigin
 
        #use the Transformation to transform the corners of the box
        v = [transform(Vector(i)) for i in v]
 
        bot =  [v[0], v[1], v[2], v[4]] # O, x, xy, y
        top =  [v[7], v[3], v[5], v[6]] # xyz, xz, z, yz
        minL = [v[0], v[4], v[6], v[5]] # O, y, yz, z
        minR = [v[0], v[5], v[3], v[1]] # O, z, xz, x
        maxL = [v[4], v[2], v[7], v[6]] # y, xy, xyz, yz
        maxR = [v[3], v[1], v[2], v[7]] # xz, x, xy, xyz
        box = [bot, minR, minL, maxR, maxL, top]
 
        return box
 
def drawBox(box, name, colourRGB):
        boxObj = []
        for side in box:
                boxObj.append(BEGIN)
                boxObj.append(LINE_STRIP)
                boxObj.append(COLOR)
                boxObj.extend(colourRGB)
                for point in side:
                        boxObj.append(VERTEX)
                        boxObj.extend(point)
                boxObj.append(END)
 
        cmd.set('auto_zoom', 0)
        cmd.load_cgo(boxObj, "%s-box" % name)
        cmd.set('auto_zoom', 1)
 
def findMinDimensionsAndRotation(coords, centroid, x, y, z):
        O = Vector(0, 0, 0)
        X = Vector(1, 0, 0)
        Y = Vector(0, 1, 0)
        Z = Vector(0, 0, 1)
 
        #Create a Transformation t = |x, y, z| . |X, Y, Z| ^ -1
        mfinal = array([array(X), array(Y), array(Z)])
        morig = array([array(x), array(y), array(z)])
        rotmatrix = matrixmultiply(morig, transpose(mfinal))
        tTotal = Transformation.Rotation(Tensor(rotmatrix))
 
        #Transform the coordinates and find the min, max dimensions
        coordArray = array([array(tTotal(Vector(c))) for c in coords])
        minDimensions = [max(coordArray[:,i]) - min(coordArray[:,i]) for i in range(3)]
 
        #Now, compose the inverse rotation used to move the bounding box to the right place
        tForward = Transformation.Rotation(Tensor(matrixmultiply(mfinal, transpose(morig))))
 
        return minDimensions, tForward
