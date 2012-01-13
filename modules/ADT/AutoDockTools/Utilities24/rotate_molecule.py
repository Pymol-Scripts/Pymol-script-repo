#!/usr/bin/env python
#$Id: rotate_molecule.py,v 1.2 2007/10/08 18:13:30 rhuey Exp $
import os 
from MolKit import Read
from MolKit.pdbWriter import PdbWriter, PdbqsWriter, PdbqWriter, PdbqtWriter
from mglutil.math.rotax import rotax
import numpy.oldnumeric as Numeric


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: rotate_molecule.py -f filename"
        print
        print "    Description of command..."
        print "        [-f]    filename"
        print "    Optional parameters:"
        print "        [-o]    alternative output filename"
        print "        (default is 'rotated_' +filename)"
        print "        [-y]    rotate around the y axis"
        print "        (default is rotation around the z axis)"
        print "        [-x]    rotate around the x axis"
        print "        (default is rotation around the z axis)"
        print "        [-u]    user-defined axis of rotation '1.0,2.0,-6.2'"
        print "        (default is rotation around the z axis)"
        print "        [-a]    angle for rotation about axis "
        print "        (default is rotation around the z axis)"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:o:xyu:a:v')

    except getopt.GetoptError, msg:
        print 'rotate_molecule.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-f: pdb_filename_stem
    filename =  None

    # optional parameters
    verbose = None
    outputfilename =  None
    rotation = 'z'
    #arbitrary axis angle for rotation
    axis = None
    angle = None

    #'f:o:v'
    for o, a in opt_list:
        print "o=", o, " a=",a
        if o in ('-f', '--f'):
            filename = a
            if verbose: print 'set filename to ', filename
            outputfilename =  'rotated_' + filename
        if o in ('-o', '--o'):
            outputfilename = a 
            if verbose: 
                print 'set output outputfilename to ', a
        if o in ('-x', '--x'):
            rotation = 'x'
            if verbose: print 'set rotation to ', rotation
        if o in ('-y', '--y'):
            rotation = 'y'
            if verbose: print 'set rotation to ', rotation
        if o in ('-u', '--u'):
            axis = a
            if verbose: print 'set user-defined axis to ', axis
        if o in ('-a', '--a'):
            angle = a
            if verbose: print 'set angle for rotation to ', angle
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not filename:
        print 'rotate_molecule: filename must be specified.'
        usage()
        sys.exit()

    mol = Read(filename)[0]
    if verbose: print 'read ', filename
    filetype = os.path.splitext(os.path.basename(filename))[1]
    if verbose: print "filetype=", filetype
    writer = None
    if filetype=='.pdbqt':
        writer = PdbqtWriter()
    elif filetype=='.pdbq':
        writer = PdbqWriter()
    elif filetype=='.pdbqs':
        writer = PdbqsWriter()
    elif filetype=='.pdb':
        writer = PdbWriter()
    else:
        print 'Sorry! Unable to write this filetype->', filetype

    center = Numeric.add.reduce(mol.allAtoms.coords)/len(mol.allAtoms)
    crds = Numeric.array(mol.allAtoms.coords)
    center = Numeric.add.reduce(crds)/len(mol.allAtoms)
    crds = crds - center
    crds = crds.tolist()
    mol.allAtoms.updateCoords(crds)
    lenCoords = len(crds)
    #rotate the atoms here
    if axis is not None and angle is not None:
        rot = (float(angle)* 3.14159/180.)%(2 * Numeric.pi)
        x = Numeric.array([0.,0.,0.])
        y = Numeric.array(map(float,axis.split(',')))
        matrix = rotax(x,y, rot)
        _ones = Numeric.ones(lenCoords, 'f')
        _ones.shape = (lenCoords,1)
        mov_coords = Numeric.concatenate((crds, _ones),1)
        newcoords = Numeric.dot(mov_coords, matrix)
        nc = newcoords[:,:3].astype('f')
        for i in range(lenCoords):
            mol.allAtoms[i]._coords[0] = nc[i].tolist()
    else:
        if rotation=='z':
            #for rotation around z-axis:
            for a in mol.allAtoms:
                a._coords[0][0] = -1.*a._coords[0][0]
                a._coords[0][1] = -1.*a._coords[0][1]
        elif rotation=='y':
            #for rotation around y-axis:
            for a in mol.allAtoms:
                a._coords[0][0] = -1.*a._coords[0][0]
                a._coords[0][2] = -1.*a._coords[0][2]
        elif rotation=='x':
            #for rotation around x-axis:
            for a in mol.allAtoms:
                a._coords[0][1] = -1.*a._coords[0][1]
                a._coords[0][2] = -1.*a._coords[0][2]
    ncrds = Numeric.array(mol.allAtoms.coords)
    ncrds = ncrds + center
    ncrds = ncrds.tolist()
    mol.allAtoms.updateCoords(ncrds)

    if writer:
        outptr = open(outputfilename, 'w')
        liglines = mol.parser.allLines
        ctr = 0
        for l in liglines:
            if l.find("ATOM")!=0 and l.find("HETATM")!=0:
                outptr.write(l)
            else:
                writer.write_atom(outptr, mol.allAtoms[ctr])
                ctr += 1
        outptr.close()


# To execute this command type:
# rotate_molecule.py -f filename [-o outputfilename -u axis -a angle to rotate] -v

