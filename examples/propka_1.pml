reinitialize
import propka
fetch 1ohr, async=0
create 1ohrB3C, 1ohr
hide everything, all
show cartoon, 1ohrB3C

cmd.wizard("mutagenesis")
cmd.do("refresh_wizard")
# To get an overview over the wizard API:
for i in dir(cmd.get_wizard()): print i

# lets mutate chain B residue 3 to CYS. (1ohrB3C)
cmd.get_wizard().set_mode("CYS")
cmd.get_wizard().do_select("/1ohrB3C//B/3")

# Select the first rotamer, which is most probable
cmd.frame(1)

# Apply the mutation
cmd.get_wizard().apply()
# Close wizard
cmd.set_wizard("done")
#OR cmd.wizard(None) 
propka resi=3
zoom /1ohrB3C//B/3

