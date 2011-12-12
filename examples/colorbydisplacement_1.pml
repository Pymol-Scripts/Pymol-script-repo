reinitialize
import colorbydisplacement

fetch 1HP1, async=0
fetch 1HPU, async=0

hide everything
### Select asymmetric units from pdb file
create O5NT, /1HP1//A
create C5NT, /1HPU//C
delete 1HP1
delete 1HPU

show cartoon, O5NT
show cartoon, C5NT

ColorByDisplacementAll O5NT, C5NT, super1=resi 26-355, super2=resi 26-355, doColor=t, doAlign=t

set_view (\
     0.094686687,   -0.390707940,    0.915631354,\
     0.809000611,   -0.505792081,   -0.299485058,\
     0.580131471,    0.769104064,    0.268191338,\
     0.000000000,    0.000000000, -280.940521240,\
    26.240486145,   46.146961212,   21.702068329,\
   231.830673218,  330.050415039,  -20.000000000 )
