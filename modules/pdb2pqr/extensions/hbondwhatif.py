"""
    Hbond extension

    Find all hydrogen bonds as determined by the DISTANCE and ANGLE cutoffs below.
    Uses PDB2PQR to determine donors and acceptors, and displays all available bonds to stdout in a WHATIF-like format.

    Authors:  Todd Dolinsky and Julie Mitchell
"""

__date__ = "17 February 2006"
__author__ = "Todd Dolinsky and Julie Mitchell"

from src.utilities import *
from src.routines import *
try:
    import Numeric
except:
    import numpy as Numeric
from math import *

ANGLE_CUTOFF = 90.0       # A - D - H(D) angle
DIST_CUTOFF = 3.30         # H(D) to A distance

def usage():
    return 'Print a list of hydrogen bonds to {output-path}.hbo'

def hbondwhatif(routines, outroot):
    """
        Print a list of hydrogen bonds.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """
    outname = outroot + ".hbo"
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
                
                if (donor.residue.chainID == acc.residue.chainID): continue
                dist = distance(donorh.getCoords(), acc.getCoords())
                if dist > DIST_CUTOFF: continue
                angle = getAngle(acc.getCoords(), donor.getCoords(), donorh.getCoords())
                if angle > ANGLE_CUTOFF: continue
				
				# the following prints the hydrogen bonds in the output format used by WHAT-IF
				# and also assigns backbone/sidechain and a psuedo h-bond score. (JC Mitchell 4/07)
				
                thisBstring='S'
                thatBstring='S'
                if (donor.isBackbone()): thisBstring='B'
                if (acc.isBackbone()): thatBstring='B'				
                if (donor.tempFactor > 60.0): continue
                if (acc.tempFactor > 60.0): continue
                score= (1.7/dist) * cos(angle * 3.142 / 180.0)
                file.write("%4d %-4s (%4d  ) %s     %-4s-> %4d %-4s (%4d  ) %s     %-4sSym=   1 Val= %6.3lf  DA=%6.2f  DHA=%6.2f (%s-%s)\n" % \
				  (donor.residue.resSeq,donor.residue.name,donor.residue.resSeq, donor.residue.chainID,donor.name,acc.residue.resSeq,acc.residue.name,acc.residue.resSeq, acc.residue.chainID,acc.name, score, dist, angle, thisBstring, thatBstring)) 

    routines.write("\n")
    file.close()

