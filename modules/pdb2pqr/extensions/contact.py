"""
    Contact extension

    Find all hydrogen bonds as determined by the DISTANCE cutoff below.
    Uses PDB2PQR to determine donors and acceptors, and displays
    all available bonds to stdout in a WHATIF-like format.

    Author:  Julie C. Mitchell
"""

__date__ = "April 2007"
__author__ = "Julie C. Mitchell"

from src.utilities import *
from src.routines import *

DIST_CUTOFF = 3.5         # max distance  

def usage():
    return 'Print a list of contacts to {output-path}.con\n'

def contact(routines, outroot):
    """
        Print a list of contacts.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """
    outname = outroot + ".con"
    file = open(outname, "w")

    # Initialize - set nearby cells, donors/acceptors
    
    cellsize = int(DIST_CUTOFF + 1.0 + 1.0) 
    protein = routines.protein
    routines.setDonorsAndAcceptors()
    routines.cells = Cells(cellsize)
    routines.cells.assignCells(protein)

    for thisatom in protein.getAtoms():

        # Grab the list of thisatoms
        if not thisatom.hdonor: continue
        thisatomhs = []
        for bond in thisatom.bonds:
            if bond.isHydrogen(): thisatomhs.append(bond)
        if thisatomhs == []: continue

      
        # For each thisatom, grab all thatatomeptors
            
        count = 0
        closeatoms = routines.cells.getNearCells(thisatom)
        for thatatom in closeatoms:
			if (thisatom.residue == thatatom.residue): continue  # comment this out to include interresidue contacts
			if (thatatom.isHydrogen()): continue
			thisdist = distance(thisatom.getCoords(), thatatom.getCoords())
			if (thisdist <= DIST_CUTOFF): 
				count = count+1
				thisBstring='S'
				thatBstring='S'
				hscore= 0.0
				if (thisatom.hdonor & thatatom.hacceptor): hscore = 1.0
				if (thisatom.hacceptor & thatatom.hdonor): hscore = 1.0
				if (thisatom.isBackbone()): thisBstring='B'
				if (thatatom.isBackbone()): thatBstring='B'
				file.write("%4d %4d %-4s (%4d  ) %s     %-4s<>%4d %-4s (%4d  ) %s     %-4s D=%6.2f  H-ene=%6.2f  Sym=  (%s-%s)\n" % \
				  (count, thisatom.residue.resSeq,thisatom.residue.name,thisatom.residue.resSeq, thisatom.residue.chainID,thisatom.name,thatatom.residue.resSeq,thatatom.residue.name,thatatom.residue.resSeq, thatatom.residue.chainID,thatatom.name, thisdist, hscore, thisBstring, thatBstring)) 
				
    routines.write("\n")
    file.close()
