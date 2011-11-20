## Author: Andreas Henschel 2006
 
from pymol import cmd
from pymol.cgo import *
 
def centerOfMass(selection):
   ## assumes equal weights (best called with "and name ca" suffix)
   model = cmd.get_model(selection)
   x,y,z=0,0,0
   for a in model.atom:
       x+= a.coord[0]
       y+= a.coord[1]
       z+= a.coord[2]
   return (x/len(model.atom), y/len(model.atom), z/len(model.atom))
 
cmd.load("/group/bioinf/Data/PDBLinks/1c7c.pdb")
cmd.select("domain", "/1c7c//A/143-283/ and name ca") ## selecting a domain
 
domainCenter=centerOfMass("domain")
 
print "Center of mass: (%.1f,%.1f,%.1f)"% domainCenter
cmd.showas("cartoon", "all")
cmd.show("spheres", "domain")
 
## Creating a sphere CGO
com = [COLOR, 1.0, 1.0, 1.0, SPHERE]+list(domainCenter) + [3.0] ## white sphere with 3A radius
cmd.load_cgo(com, "CoM")
 
cmd.zoom("1c7c", 1.0)
cmd.center("domain")
 
#ah@bioinfws19:~/Projects/PyMOL$ pymol -qc centerOfMass4.py
#Center of mass: (-1.0,24.5,48.2)
#ah@bioinfws19:~/Projects/PyMOL$
