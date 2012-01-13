reinitialize
import rotkit
 
fetch 1HP1, async=0
show_as cartoon, 1HP1
show_as sticks, 1HP1 and resn ATP
set auto_zoom, off
 
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
 
####################### Create rotate states of dye atoms ###################
##### First mutate, the mutate functions take 0.2 seconds, so we put in a refesh command to wait for everything is done
rotkit.mutate("1HP1", chain="A", resi=308, target="CYS", mutframe=1)
cmd.refresh()
rotkit.mutate("1HP1", chain="A", resi=513, target="CYS", mutframe=1)
cmd.refresh()
 
##### Create simulated dye movement atoms
pseudoatom Donor, vdw=0.5
pseudoatom Acceptor, vdw=0.5
show spheres, Donor or Acceptor 
rotkit.toline("1HP1 and resi 308 and name CA","1HP1 and resi 308 and name SG","Donor","Donor",15.0)
rotkit.toline("1HP1 and resi 513 and name CA","1HP1 and resi 513 and name SG","Acceptor","Acceptor",15.0)
 
python
Dye_ang_incr = 6
Donor_angle_range = range(0,359,Dye_ang_incr)
Acceptor_angle_range = range(0,359,Dye_ang_incr)
nrstates = len(Donor_angle_range)+1
Donor_states = 1
Acceptor_states = 1
for Donor_angle in Donor_angle_range:
    Donor_states += 1
    Donor_angle_name="Donor_%s"%(Donor_angle)
    cmd.create(Donor_angle_name,"Donor")
    rotkit.rotateline("1HP1 and resi 308 and name CA","1HP1 and resi 308 and name CB",Donor_angle,Donor_angle_name)
    # Save it as states in Donor
    cmd.create("Donor",Donor_angle_name,1,Donor_states)
    cmd.create("Donor",Donor_angle_name,1,2*nrstates-Donor_states)
    cmd.group("All_Donors",Donor_angle_name)
for Acceptor_angle in Acceptor_angle_range:
    Acceptor_states += 1
    Acceptor_angle_name="Acceptor_%s"%(Acceptor_angle)
    cmd.create(Acceptor_angle_name,"Acceptor")
    rotkit.rotateline("1HP1 and resi 513 and name CA","1HP1 and resi 513 and name CB",Acceptor_angle,Acceptor_angle_name)
    # Save it as states in Acceptor
    cmd.create("Acceptor",Acceptor_angle_name,1,Acceptor_states)
    cmd.create("Acceptor",Acceptor_angle_name,1,2*nrstates-Acceptor_states)
    cmd.group("All_Acceptors",Acceptor_angle_name)
python end
disable All_Donors
disable All_Acceptors
cmd.create("Donor","All_Donors",1,1)
cmd.create("Acceptor","All_Acceptors",1,1)
mplay
