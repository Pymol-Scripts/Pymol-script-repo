#!/usr/bin/python
 
#
# This simple script will filter through all PDBs in a directory, and for each one
# save all the ligands/heterotoms (that aren't waters) to their own file.  This
# script operates at the level of molecules, not residues, atoms, etc.  Thus, if
# you have a ligand that PyMOL is treating as ONE residue, but is actually two
# separate molecules, or a molecule and an atom, then you will get multiple files.
#
 
from glob import glob
from os import path
from pymol import stored
from string import zfill
 
theFiles = glob("../*.pdb");
 
for f in theFiles:
    # load the file
    cmd.load(f);
    # remove the protein and waters
    cmd.remove("polymer or resn HOH");
 
    cmd.select("input", "all")
    cmd.select("processed", "none")
    mol_cnt = 0
 
    while cmd.count_atoms("input"):
        # filter through the selections, updating the lists
        cmd.select("current","bymolecule first input")
        cmd.select("processed","processed or current")
        cmd.select("input","input and not current")
 
        # prepare the output parameters
        curOut = path.basename(f).split(".")[0] + "_" + zfill(str(mol_cnt),5) + "_het.pdb"
        curSel = "current"
 
        # save the file
        cmd.save( curOut, curSel );
        print "Saved " + curSel + " to " + curOut
 
        mol_cnt = mol_cnt + 1;
 
    # remove all to move to next molecule
    cmd.delete("*");        
 
