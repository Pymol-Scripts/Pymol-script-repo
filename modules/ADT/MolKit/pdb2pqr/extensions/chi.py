"""
    Chi extension

    Print the backbone chi angle for each residue in the structure.
    Chi angle is determined by the coordinates of the N, CA, CB (if
    available), and CG/OG/SG atoms (if available).

    Author:  Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Todd Dolinsky"

from src.utilities import *
from src.routines import *

def usage():
    str  = "        --chi         :  Print the per-residue backbone chi\n"
    str += "                         angle to {output-path}.chi\n"
    return str

def chi(routines, outroot):
    """
        Print the list of psi angles

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".chi"
    file = open(outname, "w")

    routines.write("\nPrinting chi angles for each residue...\n")
    routines.write("Residue     chi\n")
    routines.write("----------------\n")
    
    # Initialize some variables

    protein = routines.protein

    for residue in protein.getResidues():
        if residue.hasAtom("N"): ncoords = residue.getAtom("N").getCoords()
        else: continue

        if residue.hasAtom("CA"): cacoords = residue.getAtom("CA").getCoords()
        else: continue

        if residue.hasAtom("CB"): cbcoords = residue.getAtom("CB").getCoords()
        else: continue

        if residue.hasAtom("CG"): gcoords = residue.getAtom("CG").getCoords()
        elif residue.hasAtom("OG"): gcoords = residue.getAtom("OG").getCoords()
        elif residue.hasAtom("SG"): gcoords = residue.getAtom("SG").getCoords()
        else: continue

        chi = getDihedral(ncoords, cacoords, cbcoords, gcoords)
        routines.write("%s\t%.4f\n" % (residue, chi))
        file.write("%s\t%.4f\n" % (residue, chi))
        
    routines.write("\n")
    file.close()
