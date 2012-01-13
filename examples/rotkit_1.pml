reinitialize
import rotkit
 
fetch 1HP1, async=0
show_as cartoon, 1HP1
show_as sticks, 1HP1 and resn ATP
 
###################### Make rotation axis #################
pseudoatom axisA, vdw=1.0
pseudoatom axisB, vdw=1.0
rotkit.toline("/1HP1//A/477/C","/1HP1//A/423/CG1","axisA","axisA",20)
rotkit.toline("/1HP1//A/423/CG1","/1HP1//A/477/C","axisB","axisB",5)
show spheres, axisA or axisB 
label axisA, "axisA" 
label axisB, "axisB" 
dist rotaxis, axisA, axisB
color green, rotaxis
set dash_width, 5
set dash_gap, 0
hide label, rotaxis
 
## Create rotate states of 1HP1
create 1HP1_rot, 1HP1, 1, 1
python
ang_incr = 1
anglerange = range(2,98,ang_incr)
nrstates = len(anglerange)+1
states = 1
for angle in anglerange:
	states += 1
	rot_1HP1 = "1HP1_rot_%s"%angle
	cmd.create(rot_1HP1,"(1HP1 and resi 363-550) or (1HP1 and resn ATP)")
	rotkit.rotateline("axisA","axisB",-(angle-1),rot_1HP1)	
	cmd.create("1HP1_rot",rot_1HP1,1,states)
	cmd.create("1HP1_rot",rot_1HP1,1,2*nrstates-states)
	cmd.delete(rot_1HP1)
python end
hide cartoon, (1HP1 and resi 363-550)
hide sticks, (1HP1 and resn ATP)
mplay
