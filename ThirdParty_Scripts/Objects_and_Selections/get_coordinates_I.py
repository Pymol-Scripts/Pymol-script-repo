# This script gets a copy of the coordinates in Python,
# rotates the object about the Z axis, and then
# updates the coordinates in the original object.
 
from pymol import cmd
 
model = cmd.get_model("pept")
for a in model.atom:
   a.coord=[ -a.coord[1], a.coord[0], a.coord[2]]
 
cmd.load_model(model,"tmp")
cmd.update("pept","tmp")
cmd.delete("tmp")

