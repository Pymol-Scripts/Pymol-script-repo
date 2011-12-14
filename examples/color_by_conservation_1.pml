reinitialize
import color_by_conservation
 
# get some kinases
fetch 1opk 3dtc 3p86 2eva 3efw, async=0
 
# turn on the sequence viewer
set seq_view

# align them into the "algn" object
for x in cmd.get_names(): cmd.align(x, "3efw and c. A", object="algn")
 
# color
color_by_conservation aln=algn, as_putty=1
