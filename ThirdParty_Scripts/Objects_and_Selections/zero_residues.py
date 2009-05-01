import pymol
 
def zero_residues(sel1,offset=0):
        """
        PURPOSE: renumbers the residues so that the first one is zero, or offset           .
        USAGE: zero_residues protName    # first residue is 0
        USAGE: zero_residues protName, 5 # first residue is 5
        EXAMPLE: zero_residues *
        """
        offset = int(offset)
 
        # variable to store the offset
        stored.first = None
        # get the names of the proteins in the selection
        names = cmd.get_names(sel1)
 
        # for each name shown
        for p in names:
                # get this offset
                cmd.iterate("first %s and polymer and n. CA" % p,"stored.first=resi")
                # don't waste time if we don't have to
                if ( stored.first == offset ):
                        continue;
                # reassign the residue numbers
                cmd.alter("%s" % p, "resi=str(int(resi)-%s)" % str(int(stored.first)-offset))
                # update pymol
                cmd.sort()
 
# let pymol know about the function
cmd.extend("zero_residues", zero_residues)

