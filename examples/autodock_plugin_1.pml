reinitialize 
fetch 3IG7, async=0
#load 3ig7.pdb, 3IG7
create cdk2, 3IG7 and polymer
create EFP, 3IG7 and organic and not resn ACE
delete 3IG7
 
hide everything, all
#h_add cdk2
#h_add EFP
show_as cartoon, cdk2
show_as sticks, EFP
util.cbay EFP
 
python
# The module "prepare_ligand4.py" does not like funny names of atoms, so we have to rename them
Ligand_prop = []
cmd.iterate("EFP", "Ligand_prop.append((resi, resn, name, elem, ID))")
for resi, resn, name, elem, ID in Ligand_prop:
    print("resi %s, resn %s, name %s, elem %s, ID %s"%(resi, resn, name, elem, ID))
    cmd.alter('%s and id %s'%(resn,ID),'name=%s%s%s'%('"',elem,'"'))
# To see the change, we rewrite the list.
Ligand_prop = []
cmd.iterate("EFP", "Ligand_prop.append((resi, resn, name, elem, ID))")
for resi, resn, name, elem, ID in Ligand_prop:
    print("resi %s, resn %s, name %s, elem %s, ID %s"%(resi, resn, name, elem, ID))
python end
select flexible, byres cdk2 within 3.5 of EFP
show sticks, flexible
util.cbag flexible
disable flexible
