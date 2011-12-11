reinitialize
fetch 1ohr, async=0
import propka
import surfaceatoms
hide everything, all

### We make it in python blocks, so pymol don't speed ahead.
python
### Se version 2 of script: http://www.pymolwiki.org/index.php/FindSurfaceResidues
# When we import a module in python, the namespace is normally: module.function
resis = surfaceatoms.surfaceatoms(cutoff=10.0)
# We dont wan't to kill the server by sending hundreds of requests. So we select some few.
resis = [resis[10],resis[20],resis[30]]
for resi in resis:
	newname="1ohr%s%sC"%(resi[0],resi[1])
	cmd.create(newname,"1ohr")
	cmd.show("cartoon","1ohr%s%sC"%(resi[0],resi[1]))
	cmd.wizard("mutagenesis")
	cmd.do("refresh_wizard")
	cmd.get_wizard().set_mode("CYS")
	selection="/%s//%s/%s"%(newname,resi[0],resi[1])
	cmd.get_wizard().do_select(selection)
	cmd.frame(1)
	cmd.get_wizard().apply()
	cmd.set_wizard("done")
	# When we import a module in python, the namespace is normally: module.function
	# And we see, that propka expect resi to be in "str" format.
	# And we don't want the logtime function
	propka.propka(resi="%s"%resi[1],logtime="")
	selection="/%s//%s/%s"%(newname,resi[0],resi[1])
	cmd.select("Mutation%s%s"%(resi[0],resi[1]),"byres %s"%(selection))
	print resi
python end
cmd.disable("all")
cmd.enable("1ohr")
cmd.zoom("1ohr")
cmd.show("cartoon","1ohr")
print resis
print("Number of surface mutations: %s"%len(resis))
print("Number of residues in protein: %s"%cmd.count_atoms("1ohr and name CA"))
