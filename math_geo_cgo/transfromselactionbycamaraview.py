# transform selection coordinates by the camera view
#
# The script answers this:
#   Thanks!
#   But translate[x,y,z] only translate the molecule.
#   What I want  is to put longest length of molecule in the X axes, the 
#   second Y axes, the third z axes.
#   Just like what orient command does which change the view of camera but 
#   not the coordinates.
#   Now I want the coordinates also change after orient it.
#
cv=list(cmd.get_view())
 
cmd.transform_selection("all", \
  cv[0:3]+[0.0]+ \
  cv[3:6]+[0.0]+ \
  cv[6:9]+[0.0]+ \
  cv[12:15]+[1.0])
 
cmd.reset()
