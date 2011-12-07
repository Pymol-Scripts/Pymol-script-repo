"""
    Saltbridge extension

    Find all salt bridges as determined by the cutoff distance below.
    Uses PDB2PQR to determine atom identities and distances, and write
    out all located salt bridges to stdout.
    
    NOTE:  A bond may be labeled BOTH hbond and salt-bridge if you use both
           options in one pdb2pqr run.  Look out for double counting.

    NOTE:  This extension currently does not support salt bridges with chain termini.

    Author:  Mike Bradley (heavily copied from Todd Dolinsky's hbond extension)
"""

__date__ = "25 August 2006"
__author__ = "Mike Bradley"

from src.utilities import *
from src.routines import *

DIST_CUTOFF = 4.0         # maximum cation to anion atom distance in angstroms

def usage():
    return 'Print a list of salt bridges to {output-path}.salt'

def salt(routines, outroot):
    """
        Print a list of salt bridges.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """
    outname = outroot + ".salt"
    file = open(outname, "w")

    routines.write("Printing salt bridge list...\n")

    # Define the potential salt bridge atoms here (the current lists are for the AMBER
    # forcefield and are not necessarily exhaustive).
    posresList = ["LYS","ARG","HIP"]
    negresList = ["GLU","ASP","CYM"]
    posatomList = ["NE","NH1","NH2","NZ","ND1","NE2",]
    negatomList = ["SG","OE1","OE2","OD1","OD2"]

    # Initialize - set nearby cells
    # The cell size adds one for the salt bridge distance, and rounds up
    
    cellsize = int(DIST_CUTOFF + 1.0 + 1.0) 
    protein = routines.protein
    routines.cells = Cells(cellsize)
    routines.cells.assignCells(protein)

    # Loop over all the atoms
    for cation in protein.getAtoms():
        # check that we've found a cation
        if cation.residue.name == "NMET":
            print "YES NMET"
        if cation.residue.name not in posresList: continue
        else:
            if cation.name not in posatomList: continue
        # For each cation, grab all potential anions in nearby cells
        closeatoms = routines.cells.getNearCells(cation)
        for anion in closeatoms:
            if cation.residue.name == anion.residue.name: continue
            if anion.residue.name not in negresList: continue
            else:
                if anion.name not in negatomList: continue
            # Do distance check
            dist = distance(cation.getCoords(), anion.getCoords())
            if dist > DIST_CUTOFF: continue
            #routines.write("Cation: %s %s\tAnion: %s %s\tsaltdist: %.2f\n" % \
            #          (cation.residue, cation.name, anion.residue, anion.name, dist)) 
            file.write("Cation: %s %s\tAnion: %s %s\tsaltdist: %.2f\n" % \
                      (cation.residue, cation.name, anion.residue, anion.name, dist))
    #routines.write("\n")
    file.close()
