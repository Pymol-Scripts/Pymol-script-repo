#!/usr/bin/env python
#$Id: write_random_state_ligand.py,v 1.3 2008/03/04 17:58:57 rhuey Exp $
import os
import random
from MolKit import Read
from MolKit.molecule import BondSet
from mglutil.math.statetocoords import StateToCoords
from AutoDockTools.Conformation import Conformation
#from MolKit.pdbWriter import PdbqWriter, PdbqtWriter




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: write_random_state_ligand.py -l filename"
        print "    Description of command..."
        print "        [-l]    ligand filename"
        print "    Optional parameters:"
        print "        [-o]    alternative output filename"
        print "        (default is filename_stem+'_random_.pdbqt')"
        print "        [-t]    translation scale"
        print "                 (default is 10)"
        print "        [-x]    axis scale"
        print "                 (default is 2)"
        print "        [-a]    angle scale"
        print "                 (default is 360)"
        print "        [-d]    torsion scale"
        print "                 (default is 90)"
        print "        [-n]    number of tries"
        print "                 (default is 20)"
        print "        [-v]    verbose output"
        print

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:o:t:x:a:d:n:vh')

    except getopt.GetoptError, msg:
        print 'write_random_state_ligand.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: pdbq(t)_filename
    filename =  None
    tscale = 10.  #translation scale
    xscale = 2.   #axis scale
    ascale = 360.   #angle scale
    dscale = 90.  #torsion scale
    ntries = 20

    # optional parameters
    verbose = None
    outputfilename =  None

    #'l:o:t:x:a:d:n:vh'
    for o, a in opt_list:
        if o in ('-l', '--l'):
            filename = a
            if verbose: print 'set filename to ', filename
            stem, ext = os.path.splitext(filename)
            stem =stem.split('_')[0]
            outputfilename =  stem + '_random' + ext
            if verbose: print " default outputfilename is ", outputfilename
        if o in ('-t', '--t'):
            tscale = float(a)
            if verbose: print 'set translation scale factor to ', tscale
        if o in ('-x', '--x'):
            xscale = a
            if verbose: print 'set axis scale factor  to ', xscale
        if o in ('-a', '--a'):
            ascale = float(a)
            if verbose: print 'set angle scale factor to ', ascale
        if o in ('-d', '--d'):
            dscale = float(a)
            if verbose: print 'set torsions scale factor to ', dscale
        if o in ('-n', '--n'):
            ntries = int(a)
            if verbose: print 'set number of tries to ', ntries
        if o in ('-o', '--o'):
            outputfilename = a 
            if verbose: print 'set output outputfilename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not filename:
        print 'write_random_state_ligand: filename must be specified.'
        usage()
        sys.exit()

    m = Read(filename)[0]
    if verbose: print 'read ', filename
    #validate specified ligand file
    assert hasattr(m, 'torTree'), "specified ligand does not have a torsion tree"
    ndihe = m.parser.keys.count('BRANCH')
    if verbose: print m.name, ' has ', ndihe, ' torsions'

    #1. prepare molecule 
    m.buildBondsByDistance()
    orig_coords = m.allAtoms.coords[:]
    m.allAtoms.addConformation(m.allAtoms.coords)
    coord_index = 1
    origin = m.getCenter()
    m.stoc = StateToCoords(m, origin, 1)
    #build reference dictionary of original bonds
    orig = len(m.allAtoms.bonds[0])
    orig_d = {}
    for a in m.allAtoms: orig_d[a]=set(a.bonds.getAtoms())

    #try a new random state up to ntries times
    # convert trans to space centered on m
    for new_try in range(ntries):
        #reset coords in working slot to original coords for new try
        m.allAtoms.updateCoords(orig_coords, ind=coord_index)
        TRANS = []
        for ind in range(3):
            #do not subtract origin here, i don't understand why not
            #TRANS.append((random.uniform(-1,1)*tscale)-origin[i])
            TRANS.append(random.uniform(-1,1)*tscale)
        if verbose: print 'translation=', TRANS
        #rot = [1,0,0,180]
        #for the rotation axis, use xscale
        rot = []
        for ind in range(3):
            rot.append(random.uniform(-1,1)*xscale)
        #for the rotation angle
        rot.append(random.uniform(-1,1)*ascale)
        if verbose: print 'axis=', rot[:3], ' and angle=', rot[-1]
        #torsions:
        dihe = []
        for ind in range(ndihe):
            dihe.append(dscale * random.uniform(-1,1))
        if verbose: print 'dihe=', dihe
        conf = Conformation(m, origin, TRANS, rot, dihe)
        new_coords = conf.getCoords()
        #verify that no new bonds would be formed for this conformation
        m.allAtoms.updateCoords(new_coords, ind=coord_index)
        #remove all original bonds and set hasBonds to 0 on all levels
        del(m.allAtoms.bonds)
        m.hasBonds=0
        for c in m.chains: c.hasBonds=0
        for r in m.chains.residues: r.hasBonds=0
        for a in m.allAtoms: 
            a.bonds = BondSet([])
            a.hasBonds = 0
        m.buildBondsByDistance()
        newLen = len(m.allAtoms.bonds[0])
        if verbose: print "originally %d bonds form; after transformation %d bonds form" %(orig, newLen)
        new_d = {}
        for a in m.allAtoms: new_d[a]=set(a.bonds.getAtoms())
        ok = True
        for a in m.allAtoms:
            if orig_d[a]!=new_d[a]:
                ok = False
                if verbose: 
                    print "bonds differ for %s: %d v %d" %(a.full_name(), len(orig_d[a]), len(new_d[a]))
        #for a in m.allAtoms:
        #    assert orig_d[a]==new_d[a], '%s bonds differ: %d v %d' %(a.full_name(), len(orig_d[a]), len(new_d[a]))
        if verbose: print "end of try ", new_try+1
        if ok: 
            print "new conformation found in",
            if new_try==0:
                print "1 try" 
            else:
                print "%d tries" %(new_try+1)
            break
    if ok:
        if verbose: 
            print 'no new contacts formed so writing file after ', new_try+1,
            if new_try==0:
                print ' try...'
            else:
                print ' tries...'
        m.parser.write_with_new_coords(new_coords, outputfilename)
        if verbose: print "wrote new conformation in ", outputfilename
    else:
        print "FAILED! no new conformation found in %d attempts!!!" %ntries


# To execute this command type:
# write_random_state_ligand.py -l ligand filename -t translation -a axisangle -d torsions [-o outputfilename] -v verbose output
