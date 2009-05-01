from pymol import cmd
from pymol import stored
 
stored.xyz = []
cmd.iterate_state(1,"pept","stored.xyz.append([x,y,z])")
 
# at this point, stored.xyz is a native Python array holding
# the coordinates, which you can modify as required
 
stored.xyz = map(lambda v:[-v[1],v[0],v[2]],stored.xyz)
 
# and now you can update the internal coordinate sets
 
cmd.alter_state(1,"pept","(x,y,z)=stored.xyz.pop(0)")

