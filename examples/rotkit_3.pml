reinitialize
#You need to make sure you are in the right dir, since we are going to make some datafiles
cd /home/tlinnet/test
import rotkit
 
fetch 1HP1, async=0
show_as cartoon, 1HP1
show_as sticks, 1HP1 and resn ATP
set auto_zoom, off
 
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
Donor_names = []
Acceptor_names = []
for Donor_angle in Donor_angle_range:
    Donor_angle_name="Donor_%s"%(Donor_angle)
    Donor_names.append([Donor_angle,Donor_angle_name])
    cmd.create(Donor_angle_name,"Donor")
    rotkit.rotateline("1HP1 and resi 308 and name CA","1HP1 and resi 308 and name CB",Donor_angle,Donor_angle_name)
    cmd.group("All_Donors",Donor_angle_name)
for Acceptor_angle in Acceptor_angle_range:
    Acceptor_angle_name="Acceptor_%s"%(Acceptor_angle)
    Acceptor_names.append([Acceptor_angle,Acceptor_angle_name])
    cmd.create(Acceptor_angle_name,"Acceptor")
    rotkit.rotateline("1HP1 and resi 513 and name CA","1HP1 and resi 513 and name CB",Acceptor_angle,Acceptor_angle_name)
    cmd.group("All_Acceptors",Acceptor_angle_name)
python end
disable All_Donors
disable All_Acceptors
cmd.create("Donor","All_Donors")
cmd.create("Acceptor","All_Acceptors")
cmd.refresh()
 
# Make a distribution for the Open case
Don_Acc_distribution = []
python
for Don in Donor_names:
    for Acc in Acceptor_names:
        distname = "%s_%s"%(Don[1],Acc[1])
        distance = cmd.dist(distname,Don[1],Acc[1])
        Don_Acc_distribution.append([Don[0], Acc[1], distance])
        cmd.delete(distname)
python end
Newdir=rotkit.createdirs("results_rotkit")
os.chdir(Newdir) 
rotkit.makehistogram(Don_Acc_distribution,dataname="Don_Acc_Open",datalistindex=2,nrbins=100,binrange=[0,0])
 
# Make a distribution for angle range
cmd.create("Acceptor_rot","All_Acceptors")
python
ang_incr = 1
anglerange = range(2,98,ang_incr)
nrstates = len(anglerange)+1
states = 1
for angle in anglerange:
    states += 1
    rot_Acceptor = "Acceptor_rot_%s"%angle
    cmd.create(rot_Acceptor,"Acceptor_rot")
    rotkit.rotateline("/1HP1//A/423/CG1","/1HP1//A/477/C",-(angle-1),rot_Acceptor)	
    cmd.create("Acceptor_rot",rot_Acceptor,1,states)
    cmd.create("Acceptor_rot",rot_Acceptor,1,2*nrstates-states)
    cmd.delete(rot_Acceptor)
    for Acc in Acceptor_names:
        rotkit.rotateline("/1HP1//A/423/CG1","/1HP1//A/477/C",(angle-1),Acc[1])
python end
