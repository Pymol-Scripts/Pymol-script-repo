reinitialize
import ccp4_contact

fetch 2c7r, async=0
remove solvent
show_as cartoon, 2c7r

python
if 'PYMOL_GIT_MOD' in os.environ:
    example_dir = os.path.join(os.path.split(os.environ['PYMOL_GIT_MOD'])[0],"files_for_examples")
    contactfile = os.path.join(example_dir,"2c7r.contact")
else:
    contactfile = "2c7r.contact"
python end

select ligands, organic
select prot, chain A
select ssDNAa, chain C
select ssDNAb, chain D
select dsDNA, chain C+D

ccp4_contact.ccp4_contact(contactfile, selName1="prot", selName2="dsDNA")

# See here to represent nuc acids
#http://www.pymolwiki.org/index.php/Examples_of_nucleic_acid_cartoons
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
color slate, dsDNA and elem C

show sticks, prot_res and prot
color raspberry, prot_res and prot
show dots,  prot_atom and prot
orient prot_res
ray
