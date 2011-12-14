reinitialize
import ccp4_pisa

fetch 2c7r, async=0
preset.pretty("2c7r")

python
if 'PYMOL_GIT_MOD' in os.environ:
    example_dir = os.path.join(os.path.split(os.environ['PYMOL_GIT_MOD'])[0],"files_for_examples")
    pisafile = os.path.join(example_dir,"2c7r.pisa")
else:
    pisafile = "2c7r.pisa"
python end

ccp4_pisa.ccp4_pisa(pisafile)

python
for selname in cmd.get_names('selections')[1:]:
    cmd.create("O_%s"%selname,selname)
    cmd.show("spheres","O_%s"%selname)
    cmd.disable("O_%s"%selname)
    cmd.group("Selections",selname)
    cmd.group("Objects","O_%s"%selname)
python end

