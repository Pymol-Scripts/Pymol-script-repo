"""
    Hbond extension

    Find all hydrogen bonds as determined by the cutoffs below.
    Uses PDB2PQR to determine donors and acceptors, and displays
    all available bonds to stdout.  Please read the "NOTE" comments
    within the code for additional important information.

    Authors:  Todd Dolinsky, Michael J Bradley
"""

__date__ = "17 February 2006"
__author__ = "Todd Dolinsky"
# NOTE: This extension edited and updated on 05 August 2008 by Michael J Bradley
# NOTE: The definitions for hbonds used below were utilized in the following study:
# Bradley MJ, Chivers PT, Baker NA. Molecular dynamics simulation of the Escherichia coli NikR protein: Equilibrium conformational fluctuations reveal inter-domain allosteric communication pathways. Journal of Molecular Biology, 378, 1155-1173, 2008.  http://dx.doi.org/10.1016/j.jmb.2008.03.010

import sys

from src.utilities import *
from src.routines import *

# NOTE: You should try different angle cutoffs and see how they affect your results
ANGLE_CUTOFF = 30.0       # A - D - H(D) angle
#ANGLE_CUTOFF = 20.0       # A - D - H(D) angle

# NOTE: This hbond extension uses a different distance definition (not the H(D) - A distance).  You should still try different values for this distance cutoff and see how it affects your results.
DIST_CUTOFF = 3.4         # D to A distance

# NOTE: If you wish to use the previous hbond.py defintiion, comment out "DIST_CUTOFF" above, uncomment the next line below, and change "ANGLE_CUTOFF" to 20.0
#DIST_CUTOFF = 3.3         # H(D) to A distance


def usage():
    return 'Print a list of hydrogen bonds to {output-path}.hbond'

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

    sys.stderr.write("Warning: New H-bonding definition is being used- see hbond.py in extensions directory for details!\n")

    file.write("# Warning: New H-bonding definition is being used- see hbond.py in extensions directory for details!\n")

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
                
                #NOTE: If you wish to use the previous hbond.py definition, toggle the "dist" variable comments below
				#dist = distance(donorh.getCoords(), acc.getCoords())
                dist = distance(donor.getCoords(), acc.getCoords())
                if dist > DIST_CUTOFF: continue
                angle = getAngle(acc.getCoords(), donor.getCoords(), donorh.getCoords())
                if angle > ANGLE_CUTOFF: continue
                routines.write("Donor: %s %s\tAcceptor: %s %s\tdist: %.2f\tAngle: %.2f\n" % \
                      (donor.residue, donor.name, acc.residue, acc.name, dist, angle)) 
                file.write("Donor: %s %s\tAcceptor: %s %s\tdist: %.2f\tAngle: %.2f\n" % \
                      (donor.residue, donor.name, acc.residue, acc.name, dist, angle))

    routines.write("\n")
    file.close()
