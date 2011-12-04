"""
    Hbond extension

    Find all hydrogen bonds as determined by the cutoffs below.
    Uses PDB2PQR to determine donors and acceptors, and displays
    all available bonds to stdout.

    Author:  Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Todd Dolinsky"

from src.utilities import *
from src.routines import *

ANGLE_CUTOFF = 20.0       # A - D - H(D) angle
DIST_CUTOFF = 3.3         # H(D) to A distance

def usage():
    str  = "        --hbond       :  Print a list of hydrogen bonds to\n"
    str += "                         {output-path}.hbond\n"
    return str

def hbond(routines, outroot):
    """
        Print a list of hydrogen bonds.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """
    outname = outroot + ".hbond"
    file = open(outname, "w")

    routines.write("Printing hydrogen bond list...\n")

    # Initialize - set nearby cells, donors/acceptors
    # The cell size adds one for the D-H(D) bond, and rounds up
    
    cellsize = int(DIST_CUTOFF + 1.0 + 1.0) 
    protein = routines.protein
    routines.setDonorsAndAcceptors()
    routines.cells = Cells(cellsize)
    routines.cells.assignCells(protein)

    for donor in protein.getAtoms():

        # Grab the list of donors
        if not donor.hdonor: continue
        donorhs = []
        for bond in donor.bonds:
            if bond.isHydrogen(): donorhs.append(bond)
        if donorhs == []: continue

      
        # For each donor, grab all acceptors
            
        closeatoms = routines.cells.getNearCells(donor)
        for acc in closeatoms:
            if not acc.hacceptor: continue
            if donor.residue == acc.residue: continue
            for donorh in donorhs:

                # Do distance and angle checks
                
                dist = distance(donorh.getCoords(), acc.getCoords())
                if dist > DIST_CUTOFF: continue
                angle = getAngle(acc.getCoords(), donor.getCoords(), donorh.getCoords())
                if angle > ANGLE_CUTOFF: continue
                routines.write("Donor: %s %s\tAcceptor: %s %s\tHdist: %.2f\tAngle: %.2f\n" % \
                      (donor.residue, donor.name, acc.residue, acc.name, dist, angle)) 
                file.write("Donor: %s %s\tAcceptor: %s %s\tHdist: %.2f\tAngle: %.2f\n" % \
                      (donor.residue, donor.name, acc.residue, acc.name, dist, angle))

    routines.write("\n")
    file.close()
