#
#
# Show side chain sticks for aromatic residues
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
 
select aromatics,(resn phe+tyr+trp+his)
show sticks, (aromatics and (!name c+n+o))
color green,aromatics
disable aromatics
set cartoon_smooth_loops,0

