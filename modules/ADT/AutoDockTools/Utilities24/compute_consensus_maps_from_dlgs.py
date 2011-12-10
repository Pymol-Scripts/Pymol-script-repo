#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_consensus_maps_from_dlgs.py,v 1.8 2007/10/09 17:30:07 annao Exp $
#
# $Id: compute_consensus_maps_from_dlgs.py,v 1.8 2007/10/09 17:30:07 annao Exp $
#
import os, glob, numpy.oldnumeric as Numeric, math
from AutoDockTools.Docking import Docking
from PyAutoDock.AutoGrid import GridMap
from MolKit.molecule import Atom, AtomSet, Bond, Molecule
from MolKit.protein import Protein,Chain,ChainSet, Residue,ResidueSet
from MolKit.pdbWriter import PdbWriter

RT = 0.831*298.
from math import log, e
debug = False

def buildPdb(map_dict, npts, name='DlgBuilt', ctr=0, outputfile='results.pdb', 
                    scale=1.0):
    if debug: print "in buildPdb: tolerance=", tolerance
    name = 'DlgBuilt'
    mol = Protein(name=name)
    mol.curChain = Chain()
    mol.chains = ChainSet([mol.curChain])
    mol.curRes = Residue()
    mol.curChain.adopt(mol.curRes)
    mol.allAtoms = AtomSet()
    mol.curRes.atoms = mol.allAtoms
    nzpts=nypts=nxpts = npts
    #nxpts, nypts, nzpts = npts
    ctr = 0
    for ADtype, m in map_dict.items():
        if debug: 
            print "PROCESSING ", ADtype, " array:", max(m.ravel()), ':', min(m.ravel())
        vals = []
        tctr = 0  #for number of each type
        for z in range(nzpts):
            for y in range(nypts):
                for x in range(nxpts):
                    val = scale * abs(m[x,y,z])
                    vals.append(val)
                    #if abs(val)>.005:
                    if val>tolerance*scale:
                        ctr += 1
                        name = ADtype + str(ctr)
                        #version3:
                        #info_lo = (xcen - numxcells*spacing,
                        #    ycen - numycells*spacing, 
                        #    zcen - numzcells *spacing)
                        #using lower back pt of cube, i think
                        #xcoord = (x-info_lo[0])/spacing
                        #ycoord = (y-info_lo[1])/spacing
                        #zcoord = (z-info_lo[2])/spacing
                        #version2:
                        xcoord = (x-numxcells)*spacing + xcen
                        ycoord = (y-numycells)*spacing + ycen
                        zcoord = (z-numzcells)*spacing + zcen
                        coords = (xcoord,ycoord,zcoord)
                        tctr += 1
                    #    #print "addAtom: name=",name,"ADtype=", ADtype," val=", val, "coords=", coords,"ctr=", ctr
                        addAtom(mol, name, ADtype, val, coords, ctr)
        print "added ",tctr, '<-', ADtype, " atoms"
        if debug:
            print ADtype, ':', tctr , ' ', ctr
    print "total atoms=", ctr
    writer = PdbWriter()
    writer.write(outputfile, mol.allAtoms, records=['ATOM'])


def addAtom(mol, name, ADtype, value, coords, ctr):
    #if debug: print "in addAtom", value,
    res = mol.chains.residues[0]
    chemicalElement = ADtype[0] #???
    childIndex = ctr - 1
    top = mol
    newAt = Atom(name=name, parent=res, top=mol, 
                    chemicalElement=chemicalElement, 
                    childIndex=childIndex)
    newAt.temperatureFactor = value
    newAt.occupancy = value
    newAt.number = ctr
    newAt.conformation = 0
    newAt._coords = [list(coords)]
    newAt.hetatm = 0
    #if debug: print "added ", name, ctr, ':', newAt.full_name(),'-', newAt.parent.children.index(newAt)
    #update allAtoms attribute of this molecule
    mol.allAtoms = mol.chains.residues.atoms


def write_grid_map(filename, spacing, npts, center, score_array,maxval=10000 ):
    if debug:
        print "in write_grid_map ", filename
    stem = os.path.basename(filename).split('.')[0]
    # open and write the file
    fptr = open(filename, 'w')
    # line 1:
    ostr = "GRID_PARAMETER_FILE " + stem + ".gpf\n"
    fptr.write(ostr)
    # line 2:
    ostr = "GRID_DATA_FILE " + stem + ".maps.fld\n"
    fptr.write(ostr)
    # line 3:
    ostr = "MACROMOLECULE " + stem + ".pdbqt\n"
    fptr.write(ostr)
    # line 4:
    ostr = "SPACING " + str(spacing) + "\n"
    fptr.write(ostr)
    # line 5:
    ostr = "NELEMENTS %d %d %d\n" % (npts[0]-1, npts[1]-1, npts[2]-1)
    fptr.write(ostr)
    # line 6:
    ostr = "CENTER %f %f %f\n" % tuple(center)
    fptr.write(ostr)
    # now write the values:
    for z in range(npts[2]):
        for y in range(npts[1]):
            for x in range(npts[0]):
                value = score_array[x,y,z]
                ostr = "%.3f\n" % (value)
                fptr.write(ostr)
    # all done...
    fptr.close()


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: compute_consensus_maps_from_dlgs.py "
        print
        print "    Description of command..."
        print "    Optional parameters:"
        print "        [-n]    number of pts in each dimension (default is 101)"
        print "        [-s]    spacing between pts (default is 1.0 )"
        print "        [-t]    tolerance (default is .005 )"
        print "        [-o]    output pdb filename"
        print "                (default is 'results.pdb')"
        print "        [-a]    use all conformations  "
        print "                (default is to use only the one with the lowest energy)"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'n:s:t:o:avh')
    except getopt.GetoptError, msg:
        print 'compute_consensus_maps_from_dlgs.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    # optional parameters
    #-n: number of pts in each dimension 
    num_pts =  101
    #-s: spacing between pts 
    spacing =  1.0
    #-t: energy cutoff for including pt in calculation 
    tolerance =  0.005
    #-o outputfilename
    #-o outputfilename
    outputfilename = "results.pdb"
    #-a  use_all_conformations
    use_all_conformations = True
    #-verbose: chatty output
    verbose = None

    #'n:s:o:avh
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-n', '--n'):
            num_pts = a
            if verbose: print 'set number of pts in each dimension to ', a
        if o in ('-s', '--s'):
            spacing = float(a)
            if verbose: print 'set spacing to ', a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-a', '--a'):
            use_all_conformations = True
            if verbose: print 'set use_all_conformations to True'
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    #read all the docking logs in current directory, one by one
    dlg_list = glob.glob('./*.dlg')
    dockings = []
    #build a list of all atom types in all dlgs
    #it is assumed that all the dockings used the same grids
    ctr = 0
    at_types = {}
    for dlg in dlg_list:
        d = Docking()
        d.readDlg(dlg)
        ctr+= 1
        print ctr, ": read ", dlg
        dockings.append(d)
        for a in d.ligMol.allAtoms:
            at_types[a.autodock_element] = 0
    if debug: print 'at_types=', at_types.keys()
    d = dockings[0]  #get grid info from the first docking
    xcen, ycen, zcen = d.dlo_list[0].parser.center_pt
    #for the output maps...
    # nxgrid=nygrid=nzgrid=npts  ??is this required??
    nxpts = nypts = nzpts = int(num_pts)/2 * 2 + 1 #ensure an odd integer
    npts = (nxpts, nypts, nzpts)
    macroStem = d.dlo_list[0].macroStem
    # build list of all atom types
    all_types  = at_types.keys()
    # setup a map for each atom type
    maps = {}
    norm_maps = {}
    for t in all_types:
        maps[t] = Numeric.zeros((nxpts, nypts, nzpts)).astype('f')
        norm_maps[t] = Numeric.zeros((nxpts, nypts, nzpts)).astype('f')

    totalEnergy = 0
    #compute the totalEnergy
    if debug: print "computing the totalEnergy"
    N = 0
    for d in dockings:
        len_clusts = len(d.clusterer.clustering_dict.keys())
        if len_clusts:
            if use_all_conformations is True:
                confs = d.ch.conformations
                N += len(confs)
            else:
                key = d.clusterer.clustering_dict.keys()[0]
                confs = [d.clusterer.clustering_dict[key][0][0]]
                N += 1
        else:
            confs = [d.ch.conformations[0]]
        for c in confs:
            d.ch.set_conformation(c)
            deltaG = c.binding_energy
            c.ddG = e**(-deltaG/RT)
            totalEnergy += c.ddG
    #compute the probability of the individual conf._pi = c.ddG/totalEnergy
    print "computing the individual conf._pi"
    for d in dockings:
        len_clusts = len(d.clusterer.clustering_dict.keys())
        #if there is a clustering
        if len_clusts:
            if use_all_conformations is True:
                confs = d.ch.conformations
            else:
                key = d.clusterer.clustering_dict.keys()[0]
                confs = [d.clusterer.clustering_dict[key][0][0]]
        else:
            #if there was no clustering
            if use_all_conformations is True:
                confs = [d.ch.conformations]
            else:
                confs = [d.ch.conformations[0]]
        for c in confs:
            c._pi = c.ddG/totalEnergy
    #compute the probability of finding an atom at a specific pt in the maps
    print "totalEnergy=", totalEnergy
    print "computing the individual atom probabilities"
    for d in dockings:
        len_clusts = len(d.clusterer.clustering_dict.keys())
        if len_clusts:
            if use_all_conformations is True:
                confs = d.ch.conformations
            else:
                key = d.clusterer.clustering_dict.keys()[0]
                confs = [d.clusterer.clustering_dict[key][0][0]]
        else:
            if use_all_conformations is True:
                confs = [d.ch.conformations]
            else:
                confs = [d.ch.conformations[0]]
        for c in confs:
            #update the coordinates
            d.ch.set_conformation(c)
            for atom in d.ligMol.allAtoms:
                #get proper atomic grid
                m = maps[atom.autodock_element]
                #go to nearest pt in grid
                #FIX #1
                x,y,z = atom.coords
                numxcells = nxpts/2  #integer division
                numycells = nypts/2  #integer division
                numzcells = nzpts/2  #integer division
                #6/22: version3
                info_lo = (xcen - numxcells*spacing,
                            ycen - numycells*spacing, 
                            zcen - numzcells *spacing)
                #using lower back pt of cube, i think
                thispt = (int((x-info_lo[0])/spacing),
                          int((y-info_lo[1])/spacing),
                          int((z-info_lo[2])/spacing))
                #4/8:version 2
                #thispt = [  int(x/spacing-xcen+ spacing/2.0 +numxcells), 
                #            int(y/spacing-ycen+ spacing/2.0 +numycells), 
                #            int(z/spacing-zcen+ spacing/2.0 +numzcells)]
                #add this conf's probability to thispt ie: c._pi
                m[thispt] = m[thispt] + c._pi
    #calculate energies based on a Boltzman distribution (~normalizing the maps)
    maxval = -RT * log(tolerance)  #based on tolerance 0.005 cutoff
    if debug: print "tolerance ", tolerance, " yields maxval=", maxval
    minval = 10000
    RTtolerance = -RT*log(tolerance)
    for t in all_types:
        if debug: print "computing ", t, " normalized map"
        for z in range(nzpts):
            for y in range(nypts):
                for x in range(nxpts):
                    value = maps[t][x,y,z]
                    if value < tolerance: #-RTlog(0.005)=1312.06
                        Er = RTtolerance
                        #Er = 0.
                        #.00001
                        #Er = 2851.
                    else:
                        Er = min(-RT*log(value), 10000)
                        if Er<minval:
                            minval = Er
                    norm_maps[t][x,y,z] = Er
    #hd map: max=1290.887.0647.. min=0.0
    #maxval overal=1312.0647... minval=565.784
    # now adjust the values to the typical range of autogrid maps: -2<->2000:
    #maxminval = max(Numeric.array(minvals)) #560
    #minminval = min(Numeric.array(minvals)) #483
    rng = int(maxval-minval)  #-483
    #rng = int(half_maxval-minminval)  #-483
    if debug: print "maxval=", maxval, '- minval=', minval, ' so rng=', rng
    #print "half_maxval=", half_maxval, '- minminmval=', minminval, ' so rng=', rng
    #TRY1: rng = int(half_maxval-minminval)
    MAPS = {}
    zero_types = []
    for t in all_types:
        #do not process zero-maps
        nmap= norm_maps[t][:]
        #if max(nmap.ravel())!=min(nmap.ravel()):
        if max(nmap.ravel())==min(nmap.ravel()):
            del(norm_maps[t])
            zero_types.append(t) 
            print t, ' is a map of all zeros'
        else:
            nmap = nmap - maxval
            for z in range(nzpts):
                for y in range(nypts):
                    for x in range(nxpts):
                        value = nmap[x,y,z]
                        #if value<10:
                        nmap[x,y,z] = value/rng
                        #nmap[x,y,z] = 2.*value/rng
            MAPS[t] = nmap
    #write out the maps
    center = d.dlo_list[0].parser.center_pt
    if debug: print "writing the grid maps"
    #for atom_type, score_array in norm_maps.items():
    for atom_type, score_array in MAPS.items():
        filename = "%s.%s_combined.map" %(macroStem, atom_type)
        write_grid_map(filename, spacing, npts, center, score_array )
    #if debug: print "building the molecule"
    ####buildPdb(norm_maps, nxpts, name='DlgBuilt', ctr=0, outputfile='results.pdb')
    buildPdb(MAPS, nxpts, name='DlgBuilt', ctr=0, outputfile='results.pdb')

# To execute this command type:
# compute_consensus_maps_from_dlgs.py 
