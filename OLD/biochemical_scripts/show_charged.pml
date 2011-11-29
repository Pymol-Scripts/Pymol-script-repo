#
# Show side chain sticks for charged residues
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
 
select pos,(resn arg+lys+his)
show sticks, (pos and !name c+n+o)
color marine,pos
disable pos
select neg,(resn glu+asp)
show sticks, (neg and !name c+n+o)
color red,neg
disable neg
set cartoon_smooth_loops,0
