reinitialize
import displacementmap

fetch 1HP1, async=0
fetch 1HPU, async=0

hide everything
### Select asymmetric units from pdb file
create O5NT, /1HP1//A
create C5NT, /1HPU//C
delete 1HP1
delete 1HPU

show cartoon, all
super O5NT and resi 26-355, C5NT and resi 26-355

### Color
set_color goldenrod1, [1.000, 0.757, 0.145]
color goldenrod1, resi 26-355
set_color darkolivegreen1, [0.792, 1.000, 0.439]
color darkolivegreen1, O5NT and resi 356-362
set_color darkolivegreen4, [0.431, 0.545, 0.239]
color darkolivegreen4, C5NT and resi 356-362
set_color chocolate3, [0.804, 0.400, 0.114]
color chocolate3, O5NT and resi 363-550
set_color purple4, [0.333, 0.102, 0.545]
color purple4, C5NT and resi 363-550

# Select active site
create active_site, resi 117+120+252+116+217+84+41+43+254
show sticks, active_site
disable active_site

# Make Cys-Cys bonds
create SS, (cys/ca+cb+sg) and byres (cys/sg and bound_to cys/sg)
show sticks, SS
color yellow, SS

# Mark water molecules
create waters, resn HOH
show nb_spheres, waters
color blue, waters
disable waters

### Make sharper, and transparent
set fog=0
set cartoon_transparency, 0.7

### Load the function
dispmap
#dispmap O5NT, C5NT, 40.0, 15.0, resi1=206, resi2=1-512.515-550
## Look for serine
#dispmap mindist=40.0, mindelta=15.0, resi1=206, resi2=330.347.350.405.412.419.457.467.533.534.539.548.336.367.383.397.439.448.490.495.501.518
#dispmap resi1=308, resi2=513
zoom O5NT
