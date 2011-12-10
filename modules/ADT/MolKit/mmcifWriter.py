
# $Header: /opt/cvs/python/packages/share1.5/MolKit/mmcifWriter.py,v 1.9 2006/04/25 22:05:20 sargis Exp $
#
# $Id: mmcifWriter.py,v 1.9 2006/04/25 22:05:20 sargis Exp $
#

from os.path import splitext, basename
from string import split, strip, digits, lower, find
from MolKit.moleculeWriter import MoleculeWriter
from MolKit.pdbWriter import PdbWriter
from MolKit.tree import TreeNode, TreeNodeSet
#from MolKit import __MGLTOOLSVersion__
import os
import sys

class MMCIFWriter(MoleculeWriter):
    """Class to write data records from a molecule tree to a cif file."""
    def write(self, filename, nodes):
        assert isinstance(nodes, TreeNode) or isinstance(nodes, TreeNodeSet)
        if isinstance(nodes, TreeNode): mol = nodes.top
        elif isinstance(nodes, TreeNodeSet): mol = nodes.top.uniq()[0]
        else: return    
        
        Filename = "%s.cif"%(splitext(filename)[0])
        file_out = open(Filename, 'w')
        file_out.write("data_"+mol.name+'\n')
        file_out.write("#\n")
        
        file_out.write("_software.name                 MGLTools\n")
        #file_out.write("_software.version              " + __MGLTOOLSVersion__ + "\n")
        file_out.write("_software.contact_author       Dr. Michel Sanner\n")
        file_out.write("_software.contact_author_email mgltools@scripps.edu\n")
        file_out.write("_software.location             http://www.scripps.edu/~sanner/software\n")
        file_out.write("_software.language             Python\n")
        file_out.write("_software.compiler_version     " + sys.version.split()[0] + "\n")
        file_out.write("_software.os                   " + sys.platform + "\n")
        file_out.write("#\n_entry.id                      "+mol.name+"\n#\n")
        try:
            a,b,c = mol.cellLength
            file_out.write("_cell.length_a                 " + a + "\n")
            file_out.write("_cell.length_b                 " + b + "\n")
            file_out.write("_cell.length_c                 " + c + "\n")
            alpha,beta,gamma = mol.cellAngles
            file_out.write("_cell.angle_alpha              " + alpha + "\n")
            file_out.write("_cell.angle_beta               " + beta + "\n")
            file_out.write("_cell.angle_gamma              " + gamma + "\n")
            file_out.write("_cell.Z_PDB                    " + mol.Zvalue + "\n")
            file_out.write("_symmetry.space_group_name_H-M " + mol.spaceGroup + "\n")
        except AttributeError:
            pass
        
        file_out.write("#\n")                  
        file_out.write("loop_\n")
        file_out.write("_atom_site.group_PDB\n")     #1
        file_out.write("_atom_site.id\n")            #2
        file_out.write("_atom_site.type_symbol\n")   #3
        file_out.write("_atom_site.label_atom_id\n")  #4
        file_out.write("_atom_site.label_comp_id\n")  #5
        file_out.write("_atom_site.label_asym_id\n") #6 
        file_out.write("_atom_site.label_seq_id\n")   #7
        file_out.write("_atom_site.Cartn_x\n")       #8
        file_out.write("_atom_site.Cartn_y\n")       #9
        file_out.write("_atom_site.Cartn_z\n")       #10 
        file_out.write("_atom_site.occupancy\n")     #11    
        file_out.write("_atom_site.B_iso_or_equiv\n")#12
        
        Atoms = mol.allAtoms
        for Atom in Atoms:
            if Atom.hetatm == 0:
                file_out.write("ATOM")                              #1
            else:
                file_out.write("HETATM")
            file_out.write(" %6s"%Atom.number)                      #2
            file_out.write(" " + Atom.element)                      #3
            file_out.write(" %5s"%Atom.name)                        #4
            file_out.write(" " + Atom.parent.type)                  #5
            if Atom.parent.parent.name == ' ':
                file_out.write(" ?")                                #6
            else:
                file_out.write(" " + Atom.parent.parent.name)
            file_out.write(" %5d"%int(Atom.parent.number))          #7
            file_out.write(" %7.3f"%Atom._coords[0][0])             #8
            file_out.write(" %7.3f"%Atom._coords[0][1])             #9
            file_out.write(" %7.3f"%Atom._coords[0][2])             #10
            file_out.write(" %6.2f"%float(Atom.occupancy))          #11
            file_out.write(" %6.2f"%float(Atom.temperatureFactor))  #12
            file_out.write("\n")
        file_out.close()

if __name__ == '__main__':
    from MolKit.mmcifParser import MMCIFParser
    parser = MMCIFParser( filename='Tests/Data/1CRN.cif' )
    print "Reading molecule"
    mol = parser.parse()
    print "Done parsing"
    SS_Data  = parser.parseSSData( mol )
    print "Done parsing secondary structure"
    writer = MMCIFWriter()
    writer.write('Tests/Data/1CRN_.cif',mol)
    print "Done"
