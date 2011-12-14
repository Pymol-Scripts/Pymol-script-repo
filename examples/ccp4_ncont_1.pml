reinitialize
import ccp4_ncont

fetch 2c7r, async=0
remove solvent
show_as cartoon, 2c7r

python
if 'PYMOL_GIT_MOD' in os.environ:
    example_dir = os.path.join(os.path.split(os.environ['PYMOL_GIT_MOD'])[0],"files_for_examples")
    ncontfile = os.path.join(example_dir,"2c7r.ncont")
else:
    ncontfile = "2c7r.ncont"
python end

select ligands, organic
select prot, chain A
select ssDNAa, chain C
select ssDNAb, chain D
select dsDNA, chain C+D

ccp4_ncont.ccp4_ncont(ncontfile, selName1="prot", selName2="dsDNA")

# See here to represent nuc acids
#http://www.pymolwiki.org/index.php/Examples_of_nucleic_acid_cartoons
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
color slate, dsDNA and elem C

show sticks, prot_res
color raspberry, prot_res
show dots,  prot_atom
show dots, dsDNA_atom
orient dsDNA_res
ray
