"""
    Psi extension

    Print the phi backbone angle for each residue in the structure.
    Psi angle is determined by the coordinates of the C(i-1), N(i), CA(i), and
    C(i). atoms.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Mike Bradley, Todd Dolinsky"

from src.utilities import *
from src.routines import *

def usage():
    str  = "        --phi         :  Print the per-residue backbone phi\n"
    str += "                         angle to {output-path}.phi\n"
    return str

def phi(routines, outroot):
    """
        Print the list of phi angles

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".phi"
    file = open(outname, "w")
    
    routines.write("\nPrinting phi angles for each residue...\n")
    routines.write("Residue     Phi\n")
    routines.write("----------------\n")
    
    # Initialize some variables

    protein = routines.protein

    for residue in protein.getResidues():
        try:
            if residue.peptideC != None:
                pepcoords = residue.peptideC.getCoords()
            else: continue
        except AttributeError: # Non amino acids
            continue
        
        if residue.hasAtom("N"): ncoords = residue.getAtom("N").getCoords()
        else: continue

        if residue.hasAtom("CA"): cacoords = residue.getAtom("CA").getCoords()
        else: continue

        if residue.hasAtom("C"): ccoords = residue.getAtom("C").getCoords()
        else: continue

        phi = getDihedral(pepcoords, ncoords, cacoords, ccoords)
        routines.write("%s\t%.4f\n" % (residue, phi))
        file.write("%s\t%.4f\n" % (residue, phi))

    routines.write("\n")
    file.close()
