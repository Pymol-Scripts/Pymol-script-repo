"""
    Psi extension

    Print the psi backbone angle for each residue in the structure.
    Psi angle is determined by the coordinates of the N(i), CA(i), C(i), N(i+1)
    atoms.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Mike Bradley, Todd Dolinsky"

from src.utilities import *
from src.routines import *

def usage():
    str  = "        --psi         :  Print the per-residue backbone psi\n"
    str += "                         angle to {output-path}.psi\n"
    return str

def psi(routines, outroot):
    """
        Print the list of psi angles

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".psi"
    file = open(outname, "w")

    routines.write("\nPrinting psi angles for each residue...\n")
    routines.write("Residue     Psi\n")
    routines.write("----------------\n")
    
    # Initialize some variables

    protein = routines.protein

    for residue in protein.getResidues():
        if residue.hasAtom("N"): ncoords = residue.getAtom("N").getCoords()
        else: continue

        if residue.hasAtom("CA"): cacoords = residue.getAtom("CA").getCoords()
        else: continue

        if residue.hasAtom("C"): ccoords = residue.getAtom("C").getCoords()
        else: continue

        try:
            if residue.peptideN != None:
                pepcoords = residue.peptideN.getCoords()
            else: continue
        except AttributeError: # Non amino acids
            continue
        
        psi = getDihedral(ncoords, cacoords, ccoords, pepcoords)
        routines.write("%s\t%.4f\n" % (residue, psi))
        file.write("%s\t%.4f\n" % (residue, psi))

    routines.write("\n")
    file.close()
