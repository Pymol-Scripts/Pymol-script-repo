#
# Show side chain sticks for hydrophilic residues
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
 
select hydrophilics,(resn arg+lys+his+glu+asp+asn+gln+thr+ser+cys)
show sticks, (hydrophilics and !name c+n+o)
color green,hydrophilics
disable hydrophilics
set cartoon_smooth_loops,0
