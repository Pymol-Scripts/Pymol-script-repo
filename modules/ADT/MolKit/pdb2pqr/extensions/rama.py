"""
    Ramachandran extension

    Print both the phi and psi angles to standard out.  See the individual
    functions for more info.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Mike Bradley, Todd Dolinsky"
  
from src.utilities import *
from src.routines import *

def usage():
    str =  "        --rama        :  Print the per-residue phi and psi\n"
    str += "                         angles to {output-path}.rama for\n"
    str += "                         Ramachandran plots\n"
    return str

def rama(routines, outroot):
    """
        Print the list of phi and psi angles for use in a Ramachandran plot.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".rama"
    file = open(outname, "w")

    routines.write("\nPrinting phi and psi angles for each residue...\n")
    routines.write("Residue        Phi          Psi\n")
    routines.write("-------------------------------\n")

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
                pepncoords = residue.peptideN.getCoords()
            else: continue

            if residue.peptideC != None:
                pepccoords = residue.peptideC.getCoords()
            else: continue
        except AttributeError: # Non amino acids
            continue

        phi = getDihedral(pepccoords, ncoords, cacoords, ccoords)
        psi = getDihedral(ncoords, cacoords, ccoords, pepncoords)
        routines.write("%s\t%.4f\t%.4f\n" % (residue, phi, psi))
        file.write("%s\t%.4f\t%.4f\n" % (residue, phi, psi))

    routines.write("\n")
    file.close()
