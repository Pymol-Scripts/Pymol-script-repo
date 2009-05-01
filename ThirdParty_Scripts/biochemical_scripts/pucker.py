#fetch 1BNA.pdb
fetch 1bna
 
#Select DNA only
#Otherwise, you will get an error for water not having sugars
select DNA, not solvent
 
#Execute pucker command
pucker DNA
 
#The output should look like this
 Phase     Amp    Pucker  Residue
161.22   55.32  C2'-endo  A 1
139.52   41.67   C1'-exo  A 2
 92.82   38.31  O4'-endo  A 3
166.35   48.47  C2'-endo  A 4
128.57   46.30   C1'-exo  A 5
126.92   49.75   C1'-exo  A 6
101.30   47.32  O4'-endo  A 7
115.62   49.23   C1'-exo  A 8
140.44   46.37   C1'-exo  A 9
145.79   53.36  C2'-endo  A 10
147.47   47.04  C2'-endo  A 11
113.80   51.69   C1'-exo  A 12
 
 Phase     Amp    Pucker  Residue
153.24   43.15  C2'-endo  B 13
128.49   45.01   C1'-exo  B 14
 67.74   43.84   C4'-exo  B 15
149.33   41.01  C2'-endo  B 16
169.27   42.30  C2'-endo  B 17
147.03   42.30  C2'-endo  B 18
116.29   47.52   C1'-exo  B 19
129.62   49.92   C1'-exo  B 20
113.61   42.93   C1'-exo  B 21
156.34   50.98  C2'-endo  B 22
116.89   44.21   C1'-exo  B 23
 34.70   45.55  C3'-endo  B 24


[edit] PYMOL API

from pymol.cgo import *    # get constants
from math import *
from pymol import cmd
 
def pucker(selection):
 
	# Author: Sean Law
	# Institute: Michigan State University
	# E-mail: slaw_(at)_msu_dot_edu
 
	sel=cmd.get_model(selection)
	first=1
	old=" "
	oldchain=" "
	residue={}
	theta={}
	count=0
	for atom in sel.atom:
		new=atom.chain+" "+str(atom.resi)
		newchain=atom.chain+" "+atom.segi
		if (not (oldchain == newchain) and first):
			print " " #Blank line to separate chain output
			print "%6s  %6s  %8s  Residue" % ("Phase", "Amp", "Pucker")
		if (not(new == old) and (not first)):
			#Check that all 5 atoms exist
			if(len(residue) == 15):
				#Construct planes
				get_dihedrals(residue,theta)
				#Calculate pucker
				info = pseudo(residue,theta)
				print info+"  "+old
			else:
				print "There is no sugar in this residue"
			if (not (oldchain == newchain)):
				print " " #Blank line to separate chain output
				print "%6s  %6s  %8s  Residue" % ("Phase", "Amp", "Pucker")
			#Clear values
			residue={}
			dihedrals={}
			theta={}
			#Store new value
			store_atom(atom,residue)
		else:
			store_atom(atom,residue)
		first=0
		old=new
		oldchain=newchain
	#Final Residue
	#Calculate dihedrals for final residue
	if (len(residue) == 15):
		#Construct planes
		get_dihedrals(residue,theta)
		#Calculate pucker for final residue
		info = pseudo(residue,theta)
		print info+"  "+old
	else:
		print "There is no sugar in this residue"
	return
 
def sele_exists(sele):
	return sele in cmd.get_names("selections");
 
def pseudo(residue,theta):
	other=2*(sin(math.radians(36.0))+sin(math.radians(72.0)))
 
	#phase=atan2((theta4+theta1)-(theta3+theta0),2*theta2*(sin(math.radians(36.0))+sin(math.radians(72.0))))
	phase=atan2((theta['4']+theta['1'])-(theta['3']+theta['0']),theta['2']*other)
	amplitude=theta['2']/cos(phase)
	phase=math.degrees(phase)
	if (phase < 0):
		phase=phase+360 # 0 <= Phase < 360
	#Determine pucker
	if (phase < 36):
		pucker = "C3'-endo"
	elif (phase < 72):
		pucker = "C4'-exo"
	elif (phase <108):
		pucker = "O4'-endo"
	elif (phase < 144):
		pucker = "C1'-exo"
	elif (phase < 180):
		pucker = "C2'-endo"
	elif (phase < 216):
		pucker = "C3'-exo"
	elif (phase < 252):
		pucker = "C4'-endo"
	elif (phase < 288):
		pucker = "O4'-exo"
	elif (phase < 324):
		pucker = "C1'-endo"
	elif (phase < 360):
		pucker = "C2'-exo"
	else:
		pucker = "Phase is out of range"
	info = "%6.2f  %6.2f  %8s" % (phase, amplitude, pucker)
	return info
 
 
def store_atom(atom,residue):
	if (atom.name == "O4'" or atom.name == "O4*"):
		residue['O4*X'] = atom.coord[0]
		residue['O4*Y'] = atom.coord[1]
		residue['O4*Z'] = atom.coord[2]
	elif (atom.name == "C1'" or atom.name == "C1*"):
		residue['C1*X'] = atom.coord[0]
		residue['C1*Y'] = atom.coord[1]
		residue['C1*Z'] = atom.coord[2]
	elif (atom.name == "C2'" or atom.name == "C2*"):
		residue['C2*X'] = atom.coord[0]
		residue['C2*Y'] = atom.coord[1]
		residue['C2*Z'] = atom.coord[2]
	elif (atom.name == "C3'" or atom.name == "C3*"):
		residue['C3*X'] = atom.coord[0]
		residue['C3*Y'] = atom.coord[1]
		residue['C3*Z'] = atom.coord[2]
	elif (atom.name == "C4'" or atom.name == "C4*"):
		residue['C4*X'] = atom.coord[0]
		residue['C4*Y'] = atom.coord[1]
		residue['C4*Z'] = atom.coord[2]
	return
 
def get_dihedrals(residue,theta):
 
	C = []
	ribose = residue.keys()
	ribose.sort()
 
	shift_up(ribose,6)
	for i in range (0,12):
		C.append(residue[ribose[i]])
	theta['0']=dihedral(C)
 
	C = []
	shift_down(ribose,3)
	for i in range (0,12):
		C.append(residue[ribose[i]])
	theta['1']=dihedral(C)
 
 
	C = []
	shift_down(ribose,3)
	for i in range (0,12):
		C.append(residue[ribose[i]])
	theta['2']=dihedral(C)
 
 
	C = []
	shift_down(ribose,3)
	for i in range (0,12):
		C.append(residue[ribose[i]])
	theta['3']=dihedral(C)
 
	C = []
	shift_down(ribose,3)
	for i in range (0,12):
		C.append(residue[ribose[i]])
	theta['4']=dihedral(C)
 
	return
 
def shift_up(list,value):
	for i in range (0,value):
		list.insert(0,list.pop())
	return
 
def shift_down(list,value):
	for i in range (0,value):
		list.insert(len(list),list.pop(0))
	return
 
def dihedral(C):
 
	DX12=C[0]-C[3]
	DY12=C[1]-C[4]
	DZ12=C[2]-C[5]
 
	DX23=C[3]-C[6]
	DY23=C[4]-C[7]
	DZ23=C[5]-C[8]
 
	DX43=C[9]-C[6];
	DY43=C[10]-C[7];
	DZ43=C[11]-C[8];
 
	#Cross product to get normal
 
	PX1=DY12*DZ23-DY23*DZ12;
	PY1=DZ12*DX23-DZ23*DX12;
	PZ1=DX12*DY23-DX23*DY12;
 
	NP1=sqrt(PX1*PX1+PY1*PY1+PZ1*PZ1);
 
	PX1=PX1/NP1
	PY1=PY1/NP1
	PZ1=PZ1/NP1
 
	PX2=DY43*DZ23-DY23*DZ43;
	PY2=DZ43*DX23-DZ23*DX43;
	PZ2=DX43*DY23-DX23*DY43;
 
	NP2=sqrt(PX2*PX2+PY2*PY2+PZ2*PZ2);
 
	PX2=PX2/NP2
	PY2=PY2/NP2
	PZ2=PZ2/NP2
 
	DP12=PX1*PX2+PY1*PY2+PZ1*PZ2
 
	TS=1.0-DP12*DP12
 
	if (TS < 0):
		TS=0
	else:
		TS=sqrt(TS)
 
	ANGLE=math.pi/2.0-atan2(DP12,TS)
 
	PX3=PY1*PZ2-PY2*PZ1
	PY3=PZ1*PX2-PZ2*PX1
	PZ3=PX1*PY2-PX2*PY1
 
	DP233=PX3*DX23+PY3*DY23+PZ3*DZ23
 
	if (DP233 > 0):
		ANGLE=-ANGLE
 
	ANGLE=math.degrees(ANGLE)
 
	if (ANGLE > 180):
		ANGLE=ANGLE-360
	if (ANGLE < -180):
		ANGLE=ANGLE+360
 
	return ANGLE
 
cmd.extend("pucker",pucker)
