#
#
#
# Show side chain sticks for hydrophobic residues
#
# (script for PyMol)
# Tom Stout, 08/05/2004
#
 
mstop
dss
hide all
#zoom all
#orient
show cartoon,all
color gray,all
 
select hydrophobes,(resn ala+gly+val+ile+leu+phe+met)
show sticks, (hydrophobes and (!name c+n+o))
color orange,hydrophobes
disable hydrophobes
set cartoon_smooth_loops,0
