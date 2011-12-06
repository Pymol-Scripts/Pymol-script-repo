########################################################################
#                       AutoLigand                         
#                Author:  Rodney M. Harris                 
#                Copyright: Rodney M. Harris TSRI 2007, 2008, 2009    
#   please site:
#   Harris R, Olson A, Goodsell D.  Automated prediction
#   of ligand-binding sites in proteins.  Proteins 2008; 70:1505-1517
#########################################################################
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/AutoLigand.py,v 1.10 2009/10/01 17:12:27 sargis Exp $
#
# $Id: AutoLigand.py,v 1.10 2009/10/01 17:12:27 sargis Exp $
import sys, string
import getopt
try:
	from Numeric import zeros
except:
	from numpy.oldnumeric import zeros
from math import sqrt
import copy
def usage():
	"Print helpful usage statement to stdout."
	print "Usage:  python AutoLigand.py -r FileBaseName -p #_of_pts"
	print
	print "    Description of command..."
	print "       -r FileBaseName = just the name part from receptor map files (i.e., FileBaseName.C.map)"
	print "       -p #_of_pts = number of fill points to use (int)"
	print "              Note: can be omitted if -a option used."
	print "    Optional parameters:"
	print "       [-a #] = number of heavy atom for ligand (#_of_pts will be set to 6x atoms)"
	print "       [-x # -y # -z #] = optional x,y,z co-ords for starting fill (float)"
	print "              when starting point is input, only one fill will be run"
	print "       [-i # -j # -k #] = optional x,y,z co-ords for second point (float)"
	print "              when second point is input, the fill will connect both points"
	print "              NOTE: the connection path has not been optimized - use with discretion"
	print "       [-f #] = number of fills to generate - default is 10"
	print "       [-e] = use the extra atom types NA, N, SA, and A"
	print "              NOTE: these results can be problematic - use with discretion"
	print "       [-m] = make a movie of output fill progress"
# process command arguments
try:
	opt_list, args = getopt.getopt(sys.argv[1:], 'r:p:a:x:y:z:i:j:k:f:em')
except getopt.GetoptError, msg:
	print 'AutoLigand.py: %s' %msg
	usage()
	sys.exit(2)
# initialize parameters
FileBaseName = None
volume = 0
x_flag = False
y_flag = False
z_flag = False
i_flag = False
j_flag = False
k_flag = False
Number_fills = 10
more_atom_types = False
movie = False
# 'r:p:a:x:y:z:i:j:k:f:xm'
for o, a in opt_list:
	if o in ('-r', '--r'):
		FileBaseName = a
	if o in ('-p', '--p'):
		volume = int(a)
	if o in ('-a', '--a'):
		no_heavy_atoms = int(a)
		volume = no_heavy_atoms * 6
	if o in ('-x', '--x'):
		x_input = float(a)
		x_flag = True
	if o in ('-y', '--y'):
		y_input = float(a)
		y_flag = True
	if o in ('-z', '--z'):
		z_input = float(a)
		z_flag = True
	if o in ('-i', '--i'):
		x2_input = float(a)
		i_flag = True
	if o in ('-j', '--j'):
		y2_input = float(a)
		j_flag = True
	if o in ('-k', '--k'):
		z2_input = float(a)
		k_flag = True
	if o in ('-f', '--f'):
		Number_fills = int(a)
	if o in ('-e', '--e'):
		more_atom_types = True
	if o in ('-m', '--m'):
		movie = True
if not FileBaseName and not volume:
	usage()
	sys.exit(2)
if not x_flag:
	start_point_flag = 0
if x_flag and y_flag and z_flag and not i_flag:
	start_point_flag = 1
if x_flag and y_flag and z_flag and i_flag and j_flag and k_flag:
	start_point_flag = 2
# GUI
if movie:
	import cPickle
	output = open(FileBaseName+'_flood.pkl', 'wb')
prev_flood = []

def save_flood(flood):
	global prev_flood
	add_set = []
	for item in flood:
		if item in prev_flood:
			prev_flood.remove(item)
		else:
			add_set.append(item)
	if prev_flood or add_set:
		cPickle.dump([prev_flood, add_set], output, -1)
	prev_flood = copy.copy(flood)

every = 3 # note, this is hard-coded - if it needs changing, the StartPointModList must be changed below
# in this version of AutoLigand, only C, OA, and HD are used - the user can also add A, N, NA, SA
Aflag = 0
Cflag = 1
Hflag = 1
Nflag = 0
NAflag = 0
Oflag = 1
Sflag = 0
CFileName = FileBaseName + '.C.map'
HFileName = FileBaseName + '.HD.map'
OFileName = FileBaseName + '.OA.map'
dFileName = FileBaseName + '.d.map'
eFileName = FileBaseName + '.e.map'
try:
	Cmap = open(CFileName,'r')
except IOError:
	Cflag = 0
try:
	Hmap = open(HFileName,'r')
except IOError:
	Hflag = 0
try:
	Omap = open(OFileName,'r')
except IOError:
	Oflag = 0
emap = open(eFileName,'r')
dmap = open(dFileName,'r')
if more_atom_types:
	Aflag = 1
	Nflag = 1
	NAflag = 1
	Sflag = 1
	AFileName = FileBaseName + '.A.map'
	NFileName = FileBaseName + '.N.map'
	NAFileName = FileBaseName + '.NA.map'
	SFileName = FileBaseName + '.SA.map'
	try:
		Amap = open(AFileName,'r')
	except IOError:
		Aflag = 0
	try:
		Nmap = open(NFileName,'r')
	except IOError:
		Nflag = 0
	try:
		NAmap = open(NAFileName,'r')
	except IOError:
		NAflag = 0
	try:
		Smap = open(SFileName,'r')
	except IOError:
		Sflag = 0
for i in range(6):
	linee = emap.readline()
	lined = dmap.readline()
	if Aflag == 1: lineA = Amap.readline()
	if Cflag == 1: lineC = Cmap.readline()
	if Hflag == 1: lineH = Hmap.readline()
	if Nflag == 1: lineN = Nmap.readline()
	if NAflag == 1: lineNA = NAmap.readline()
	if Oflag == 1: lineO = Omap.readline()
	if Sflag == 1: lineS = Smap.readline()
	if i == 3:
		word = string.split(linee)
		spacing = string.atof(word[1])
	if i == 4:
		word = string.split(linee)
		nelementsx = string.atoi(word[1])
		nelementsy = string.atoi(word[2])
		nelementsz = string.atoi(word[3])
	if i == 5:
		word = string.split(linee)
		centerx = string.atof(word[1])
		centery = string.atof(word[2])
		centerz = string.atof(word[3])
Alist = []
Clist = []
Hlist = []
Nlist = []
NAlist = []
Olist = []
Slist = []
MINlist = []
MINatom = []
# charge parameters
qA = 0.01
qC = 0.0
qH = 0.163
qN = -0.24
qNA = -0.25
qO = -0.395
qS = -0.180
nx = nelementsx + 1
ny = nelementsy + 1
nz = nelementsz + 1
total_points = nx * ny * nz
i = 0
print 'total grid points = ', total_points

# GUI
if movie:
	print "output progress file = data.pkl"
	data = [int(nx/2), int(ny/2), int(nz/2),
		    centerx, centery, centerz,
		    spacing]
	cPickle.dump(data, output, -1)

while i < total_points:
	linee = emap.readline()
	worde = string.split(linee)
	eMapVal = string.atof(worde[0])
	lined = dmap.readline()
	wordd = string.split(lined)
	dMapVal = string.atof(wordd[0])
	if Aflag == 1:
		lineA = Amap.readline()
		word = string.split(lineA)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qA * eMapVal) + (abs(qA) * dMapVal)
		datum = rawdatum + (abs(qA) * dMapVal)
		Alist.append(datum)
	else: Alist.append(1000000.)
	if Cflag == 1:
		lineC = Cmap.readline()
		word = string.split(lineC)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qC * eMapVal) + (abs(qC) * dMapVal)
		datum = rawdatum + (abs(qC) * dMapVal)
		Clist.append(datum)
	else: Clist.append(900000.)
	# the Clist is set lower than the others as a default
	if Hflag == 1:
		lineH = Hmap.readline()
		word = string.split(lineH)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qH * eMapVal) + (abs(qH) * dMapVal)
		datum = rawdatum + (abs(qH) * dMapVal)
		Hlist.append(datum)
	else: Hlist.append(1000000.)
	if Nflag == 1:
		lineN = Nmap.readline()
		word = string.split(lineN)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qN * eMapVal) + (abs(qN) * dMapVal)
		datum = rawdatum + (abs(qN) * dMapVal)
		Nlist.append(datum)
	else: Nlist.append(1000000.)
	if NAflag == 1:
		lineNA = NAmap.readline()
		word = string.split(lineNA)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qNA * eMapVal) + (abs(qNA) * dMapVal)
		datum = rawdatum + (abs(qNA) * dMapVal)
		NAlist.append(datum)
	else: NAlist.append(1000000.)
	if Oflag == 1:
		lineO = Omap.readline()
		word = string.split(lineO)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qO * eMapVal) + (abs(qO) * dMapVal)
		datum = rawdatum + (abs(qO) * dMapVal)
		Olist.append(datum)
	else: Olist.append(1000000.)
	if Sflag == 1:
		lineS = Smap.readline()
		word = string.split(lineS)
		rawdatum = string.atof(word[0])
		# datum = string.atof(word[0])
		# datum = rawdatum + (qS * eMapVal) + (abs(qS) * dMapVal)
		datum = rawdatum + (abs(qS) * dMapVal)
		Slist.append(datum)
	else: Slist.append(1000000.)
	# calculate which atom is minimum value and load in list
	# atomtypes: NA = 7, S=6, A=5, O=4, N=3, C=2, H=1
	minvalue = min([Slist[i],Alist[i],Clist[i],Hlist[i],Nlist[i],NAlist[i],Olist[i]])
	if minvalue == NAlist[i]:
		atomtype = 7
	if minvalue == Slist[i]:
		atomtype = 6
	if minvalue == Alist[i]:
		atomtype = 5
	if minvalue == Nlist[i]:
		atomtype = 3
	if minvalue == Clist[i]:
		atomtype = 2
	if minvalue == Hlist[i]:
		atomtype = 1
	if minvalue == Olist[i]:
		atomtype = 4
	MINlist.append(minvalue)
	MINatom.append(atomtype)
	i = i + 1
# now load up a 3D array with zeros and fill with best affinity values
# load a second array with the best atom type
bestmap = zeros((nx,ny,nz),'float')
bestatom = zeros((nx,ny,nz), 'int')
listinc = 0
for k in range(nz):
	for j in range(ny):
		for i in range(nx):
			bestmap[i,j,k] = MINlist[listinc]
			bestatom[i,j,k] = MINatom[listinc]
			listinc = listinc + 1
# define function for finding neighbors
def nborfinder(x,nx,y,ny,z,nz,bestmap,bestatom,nbor,flood):
	"This funtion finds orthogonal neighbor points"
	if x+1 <= nx-1:
		node = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
		if node not in nbor and node not in flood:
			nbor.append(node)
	if x-1 >= 0:
		node = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
		if node not in nbor and node not in flood:
			nbor.append(node)
	if y+1 <= ny-1:
		node = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
		if node not in nbor and node not in flood:
			nbor.append(node)
	if y-1 >= 0:
		node = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
		if node not in nbor and node not in flood:
			nbor.append(node)
	if z+1 <= nz-1:
		node = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
		if node not in nbor and node not in flood:
			nbor.append(node)
	if z-1 >= 0:
		node = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
		if node not in nbor and node not in flood:
			nbor.append(node)
# define second function for finding neighbors with no steric hindrence
def nborfinder2(x,nx,y,ny,z,nz,bestmap,bestatom,nbor,flood):
	"This funtion finds orthogonal neighbor points with values <= 0"
	if x+1 <= nx-1:
		node = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
		if node[0] <= 0.0:
			if node not in nbor and node not in flood:
				nbor.append(node)
	if x-1 >= 0:
		node = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
		if node[0] <= 0.0:
			if node not in nbor and node not in flood:
				nbor.append(node)
	if y+1 <= ny-1:
		node = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
		if node[0] <= 0.0:
			if node not in nbor and node not in flood:
				nbor.append(node)
	if y-1 >= 0:
		node = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
		if node[0] <= 0.0:
			if node not in nbor and node not in flood:
				nbor.append(node)
	if z+1 <= nz-1:
		node = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
		if node[0] <= 0.0:
			if node not in nbor and node not in flood:
				nbor.append(node)
	if z-1 >= 0:
		node = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
		if node[0] <= 0.0:
			if node not in nbor and node not in flood:
				nbor.append(node)
# define function to set state for contiguous
def setstate(worstpoint,flood):
	"This function determines if a point is removable or not"
	# return a value for contig, 0 for removable and 1 for non-removable
	# create a dictionary of coords in flood and assign value to 1
	floodcoords = {}
	for fl in flood:
		floodcoords[fl[1],fl[2],fl[3]] = 1
	F = [0,0,0,0,0,0]
	tag = [0,0,0,0,0,0]
	nborcoords = []
	neighbors = 0
 	contig = 0  # defult is removable
	planes = 0
	x = worstpoint[1]
	y = worstpoint[2]
	z = worstpoint[3]
	# check for orthogonal neighbors and set F to 1 if there
	if floodcoords.has_key((x+1,y,z)):
		F[0] = 1
		neighbors = neighbors + 1
		nborcoords.append((x+1,y,z))
	if floodcoords.has_key((x-1,y,z)):
		F[3] = 1
		neighbors = neighbors + 1
		nborcoords.append((x-1,y,z))
	if floodcoords.has_key((x,y+1,z)):
		F[1] = 1
		neighbors = neighbors + 1
		nborcoords.append((x,y+1,z))
	if floodcoords.has_key((x,y-1,z)):
		F[4] = 1
		neighbors = neighbors + 1
		nborcoords.append((x,y-1,z))
	if floodcoords.has_key((x,y,z+1)):
		F[2] = 1
		neighbors = neighbors + 1
		nborcoords.append((x,y,z+1))
	if floodcoords.has_key((x,y,z-1)):
		F[5] = 1
		neighbors = neighbors + 1
		nborcoords.append((x,y,z-1))
	# check for planes, planes must tag each neighbor or point non-removable
	if floodcoords.has_key((x,y-1,z+1)) and F[4] == 1 and F[2] == 1:
		planes = planes + 1
		tag[4] = 1
		tag[2] = 1
	if floodcoords.has_key((x,y+1,z+1)) and F[1] == 1 and F[2] == 1:
		planes = planes + 1
		tag[1] = 1
		tag[2] = 1
	if floodcoords.has_key((x,y-1,z-1)) and F[4] == 1 and F[5] == 1:
		planes = planes + 1
		tag[4] = 1
		tag[5] = 1
	if floodcoords.has_key((x,y+1,z-1)) and F[1] == 1 and F[5] == 1:
		planes = planes + 1
		tag[1] = 1
		tag[5] = 1
	if floodcoords.has_key((x+1,y,z+1)) and F[2] == 1 and F[0] == 1:
		planes = planes + 1
		tag[2] = 1
		tag[0] = 1
	if floodcoords.has_key((x-1,y,z+1)) and F[2] == 1 and F[3] == 1:
		planes = planes + 1
		tag[2] = 1
		tag[3] = 1
	if floodcoords.has_key((x+1,y,z-1)) and F[0] == 1 and F[5] == 1:
		planes = planes + 1
		tag[0] = 1
		tag[5] = 1
	if floodcoords.has_key((x-1,y,z-1)) and F[3] == 1 and F[5] == 1:
		planes = planes + 1
		tag[3] = 1
		tag[5] = 1
	if floodcoords.has_key((x+1,y-1,z)) and F[0] == 1 and F[4] == 1:
		planes = planes + 1
		tag[0] = 1
		tag[4] = 1
	if floodcoords.has_key((x+1,y+1,z)) and F[0] == 1 and F[1] == 1:
		planes = planes + 1
		tag[0] = 1
		tag[1] = 1
	if floodcoords.has_key((x-1,y-1,z)) and F[3] == 1 and F[4] == 1:
		planes = planes + 1
		tag[3] = 1
		tag[4] = 1
	if floodcoords.has_key((x-1,y+1,z)) and F[1] == 1 and F[3] == 1:
		planes = planes + 1
		tag[1] = 1
		tag[3] = 1
	# now check for removable point and set state flag
	tagtotal = tag[0] + tag[1] + tag[2] + tag[3] + tag[4] + tag[5]
	if neighbors == 1:
		contig = 0  # state is set to removable
	if neighbors == 2:
		if planes >= 1 and tagtotal == 2:
			contig = 0
		else: contig = 1
	if neighbors == 3:
		if planes >= 2 and tagtotal == 3:
			contig = 0
		else: contig = 1
	if neighbors == 4:
		if planes >= 3 and tagtotal == 4:
			contig = 0
		else: contig = 1
	if neighbors == 5:
		if planes >= 4 and tagtotal == 5:
			contig = 0
		else: contig = 1
	if neighbors == 6:
		if planes >= 5 and tagtotal == 6:
			contig = 0
		else: contig = 1
	# now check if non-removable point is in a circuit and can be removed
	if contig == 1:
		# note, this can only be true if 2 or more neighbors
		# load up a test coordinate dictionary
		testfloodcoords = {}
		for tfl in flood:
			testfloodcoords[tfl[1],tfl[2],tfl[3]] = 1
		# start with the last coord added to nborcoords
		coord = nborcoords.pop()
		del testfloodcoords[(x,y,z)]  # remove starting point
		# load que with neighbors of the active point: ap
		Q = []
		apx = coord[0]
		apy = coord[1]
		apz = coord[2]
		if testfloodcoords.has_key((apx+1,apy,apz)):
			d = sqrt(((x-(apx+1))**2)+((y-apy)**2)+((z-apz)**2))
			Q.append((d,apx+1,apy,apz))
		if testfloodcoords.has_key((apx-1,apy,apz)):
			d = sqrt(((x-(apx-1))**2)+((y-apy)**2)+((z-apz)**2))
			Q.append((d,apx-1,apy,apz))
		if testfloodcoords.has_key((apx,apy+1,apz)):
			d = sqrt(((x-apx)**2)+((y-(apy+1))**2)+((z-apz)**2))
			Q.append((d,apx,apy+1,apz))
		if testfloodcoords.has_key((apx,apy-1,apz)):
			d = sqrt(((x-apx)**2)+((y-(apy-1))**2)+((z-apz)**2))
			Q.append((d,apx,apy-1,apz))
		if testfloodcoords.has_key((apx,apy,apz+1)):
			d = sqrt(((x-apx)**2)+((y-apy)**2)+((z-(apz+1))**2))
			Q.append((d,apx,apy,apz+1))
		if testfloodcoords.has_key((apx,apy,apz-1)):
			d = sqrt(((x-apx)**2)+((y-apy)**2)+((z-(apz-1))**2))
			Q.append((d,apx,apy,apz-1))
		del testfloodcoords[(apx,apy,apz)]  # delete the active point
		# now traverse the que and check for nborcoords, remove points as you go
		while 1:
			if len(Q) == 0:
				break
			Q.sort()
			qcoord = Q.pop(0)
			apx = qcoord[1]
			apy = qcoord[2]
			apz = qcoord[3]
			del testfloodcoords[(apx,apy,apz)]  # delete the new active point
			# check to see if this point is the last untouched neighbor point
			if nborcoords.count((apx,apy,apz)) == 1:
				nborcoords.remove((apx,apy,apz))
			if len(nborcoords) == 0:
				contig = 0  # all neighbor points still conected, so removable
				break
			if testfloodcoords.has_key((apx+1,apy,apz)):
				d = sqrt(((x-(apx+1))**2)+((y-apy)**2)+((z-apz)**2))
				if Q.count((d,apx+1,apy,apz)) == 0:
					Q.append((d,apx+1,apy,apz))
			if testfloodcoords.has_key((apx-1,apy,apz)):
				d = sqrt(((x-(apx-1))**2)+((y-apy)**2)+((z-apz)**2))
				if Q.count((d,apx-1,apy,apz)) == 0:
					Q.append((d,apx-1,apy,apz))
			if testfloodcoords.has_key((apx,apy+1,apz)):
				d = sqrt(((x-apx)**2)+((y-(apy+1))**2)+((z-apz)**2))
				if Q.count((d,apx,apy+1,apz)) == 0:
					Q.append((d,apx,apy+1,apz))
			if testfloodcoords.has_key((apx,apy-1,apz)):
				d = sqrt(((x-apx)**2)+((y-(apy-1))**2)+((z-apz)**2))
				if Q.count((d,apx,apy-1,apz)) == 0:
					Q.append((d,apx,apy-1,apz))
			if testfloodcoords.has_key((apx,apy,apz+1)):
				d = sqrt(((x-apx)**2)+((y-apy)**2)+((z-(apz+1))**2))
				if Q.count((d,apx,apy,apz+1)) == 0:
					Q.append((d,apx,apy,apz+1))
			if testfloodcoords.has_key((apx,apy,apz-1)):
				d = sqrt(((x-apx)**2)+((y-apy)**2)+((z-(apz-1))**2))
				if Q.count((d,apx,apy,apz-1)) == 0:
					Q.append((d,apx,apy,apz-1))
	return contig
# find the center node
xcent = int(nx/2)
ycent = int(ny/2)
zcent = int(nz/2)
if start_point_flag == 0:
	print 'done loading affinity maps, now track through all points and find (ten) best start points'
if start_point_flag == 2:
	print 'done loading affinity maps, now generate fill that contains both input points'
# first, outside mesh loop, set up a coord mods list to find 1A sphere around flood points for vol calc
CoordModsList = []
cpr = 1.0/spacing  # cpr = coord point radii
vpr = int(cpr)  # vpr = volume point radii - the number of gridpoints in radius
for vz in range(-vpr,vpr+1):
	for vy in range(-vpr,vpr+1):
		for vx in range(-vpr,vpr+1):
			voldist = sqrt((vx**2)+(vy**2)+(vz**2))
			if voldist <= cpr:
				CoordModsList.append([vx,vy,vz])
# also, make a 6 A sphere for the PP migration
# NOTE: if a 1A spacing is not used, then the PP method will search a smaller radii
CoordModsList6A = []
cpr = 6.0/spacing  # cpr = coord point radii
vpr = int(cpr)  # vpr = volume point radii - the number of gridpoints in radius
for vz in range(-vpr,vpr+1):
	for vy in range(-vpr,vpr+1):
		for vx in range(-vpr,vpr+1):
			voldist = sqrt((vx**2)+(vy**2)+(vz**2))
			if voldist <= cpr:
				CoordModsList6A.append([vx,vy,vz])
if start_point_flag == 0:
	mesh = []
	fullgrid = []
	fullmesh = []
	nonoverlapfills = []
	for k in range(0,nz,every):
		for j in range(0,ny,every):
			for i in range(0,nx,every):
				best = [i,j,k]
				fullgrid.append(best)
	# generate coords in cube about point in fullgrid
	StartPointModList = []
	for k in range(-1,2):
		for j in range(-1,2):
			for i in range(-1,2):
				StartPointModList.append([i,j,k])
	for fg in fullgrid:
		Shell = []
		for SPML in StartPointModList:
			x = fg[0]+SPML[0]
			if x >= nx or x < 0:
				continue
			y = fg[1]+SPML[1]
			if y >= ny or y < 0:
				continue
			z = fg[2]+SPML[2]
			if z >= nz or z < 0:
				continue
			node = [bestmap[x,y,z],x,y,z,bestatom[x,y,z]]
			if node[0] >= 0.0: # don't start inside receptor or far away
				continue
			if node[4] == 1:   # don't start with H atoms
				continue
			Shell.append(node)
		if len(Shell) == 0:
			continue
		# pick best point in shell
		Shell.sort()
		bestShell = Shell[0]
		spx = bestShell[1]
		spy = bestShell[2]
		spz = bestShell[3]
		# the start point should now be:  [spx,spy,spz]
		fill = []
		startpoint = [bestmap[spx,spy,spz],spx,spy,spz,bestatom[spx,spy,spz]]
		fill.append(startpoint)
		# now generate neighbors - save in growing list of all neighbors
		nbor = []
		TotalEnergy = startpoint[0]
		for i in range((volume - 1)):
			nborfinder(spx,nx,spy,ny,spz,nz,bestmap,bestatom,nbor,fill)
			# now sort neighbors list for min energy and add to fill list
			nbor.sort()
			node = nbor.pop(0)
			fill.append(node)
			TotalEnergy = TotalEnergy + node[0]
			spx = node[1]
			spy = node[2]
			spz = node[3]
		# now calculate volume of fill with 1A radii
		TotalVolume = []
		for fl in fill:
			for CML in CoordModsList:
				point = [fl[1]+CML[0],fl[2]+CML[1],fl[3]+CML[2]]
				if point not in TotalVolume:
					TotalVolume.append(point)
		NumberVolumePoints = len(TotalVolume)
		volangstrom = NumberVolumePoints * (spacing**3)
		epv = TotalEnergy/volangstrom
		fillpoint = [epv,startpoint[1],startpoint[2],startpoint[3],fill]
		fullmesh.append(fillpoint)
	fullmesh.sort()
	# now check if fills overlap and remove those that overlap by 80 percent
	# load the first fill into nonoverlapfills
	fp = fullmesh[0][4]
	for point in fp:
		nonoverlapfills.append(point)
	eightypercent = int(volume *0.8)
	tempmesh = fullmesh[1:]
	for fm in tempmesh:
		fp = fm[4]
		overlapcount = 0
		for point in fp:
			if point in nonoverlapfills:
				overlapcount = overlapcount + 1
		if overlapcount >= eightypercent:
			# remove fm from fullmesh
			fullmesh.remove(fm)
		else:
			for point in fp:
				nonoverlapfills.append(point)
	fmrange = len(fullmesh)
	if fmrange > Number_fills:
		fmrange = Number_fills
	for fm in range(fmrange):
		mesh.append(fullmesh[fm])
		testpoint = fullmesh[fm]
		x = (testpoint[1] - xcent)*spacing + centerx
		y = (testpoint[2] - ycent)*spacing + centery
		z = (testpoint[3] - zcent)*spacing + centerz
if start_point_flag == 1:
	# generate a startpoint from the input xyz
	mesh = []
	testmesh = []
	StartPointModList = []
	for k in range(-4,5):
		for j in range(-4,5):
			for i in range(-4,5):
				dist = sqrt((i**2)+(j**2)+(k**2))
				StartPointModList.append([dist,i,j,k])
	xpt = xcent + int((x_input - centerx)/spacing)
	ypt = ycent + int((y_input - centery)/spacing)
	zpt = zcent + int((z_input - centerz)/spacing)
	for SPML in StartPointModList:
		if xpt+SPML[1] >= nx-1 or xpt+SPML[1] <= 0:
			continue
		if ypt+SPML[2] >= ny-1 or ypt+SPML[2] <= 0:
			continue
		if zpt+SPML[3] >= nz-1 or zpt+SPML[3] <= 0:
			continue
		if bestmap[xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]] >= 0.0:
			continue
		if bestatom[xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]] == 1:
			continue
		fillpoint = [SPML[0],xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]]
		testmesh.append(fillpoint)
	testmesh.sort()
	if len(testmesh) == 0:
		print 'no good starting points within 4A of your chosen point'
		sys.exit(2)
	goodfillpoint = testmesh[0]
	mesh.append(goodfillpoint)
if start_point_flag == 2:
	# fill mesh with one dummy point so that the meshpoint loop only runs once
	mesh = []
	dummyfillpoint = [0.,0,0,0]
	mesh.append(dummyfillpoint)
# now start with the startpoints in mesh
itmesh = 0
unusablepointcount = 0
badHstartFlag = 0
PDBlist = []
OUTlist = []
OUTcount = 0
for meshpoint in mesh:
	itmesh = itmesh + 1
	if start_point_flag == 2:
		# load in the two fixed start points and make sure they are valid
		StartPointModList = []
		for k in range(-4,5):
			for j in range(-4,5):
				for i in range(-4,5):
					dist = sqrt((i**2)+(j**2)+(k**2))
					StartPointModList.append([dist,i,j,k])
		# generate the first startpoint from the input xyz
		testmesh = []
		xpt = xcent + int((x_input - centerx)/spacing)
		ypt = ycent + int((y_input - centery)/spacing)
		zpt = zcent + int((z_input - centerz)/spacing)
		for SPML in StartPointModList:
			if xpt+SPML[1] >= nx-1 or xpt+SPML[1] <= 0:
				continue
			if ypt+SPML[2] >= ny-1 or ypt+SPML[2] <= 0:
				continue
			if zpt+SPML[3] >= nz-1 or zpt+SPML[3] <= 0:
				continue
			if bestmap[xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]] > 0.0:
				continue
			fillpoint = [SPML[0],xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]]
			testmesh.append(fillpoint)
		testmesh.sort()
		if len(testmesh) == 0:
			print 'no good starting points within 4A of 1st input point'
			sys.exit(2)
		gpf1 = testmesh[0] # good fill point coords
		loadpoint1 = [bestmap[gpf1[1],gpf1[2],gpf1[3]],gpf1[1],gpf1[2],gpf1[3],bestatom[gpf1[1],gpf1[2],gpf1[3]]]
		print 'first user point loaded'
		# generate the second startpoint from the input xyz
		testmesh = []
		xpt = xcent + int((x2_input - centerx)/spacing)
		ypt = ycent + int((y2_input - centery)/spacing)
		zpt = zcent + int((z2_input - centerz)/spacing)
		for SPML in StartPointModList:
			if xpt+SPML[1] >= nx-1 or xpt+SPML[1] <= 0:
				continue
			if ypt+SPML[2] >= ny-1 or ypt+SPML[2] <= 0:
				continue
			if zpt+SPML[3] >= nz-1 or zpt+SPML[3] <= 0:
				continue
			if bestmap[xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]] > 0.0:
				continue
			fillpoint = [SPML[0],xpt+SPML[1],ypt+SPML[2],zpt+SPML[3]]
			testmesh.append(fillpoint)
		testmesh.sort()
		if len(testmesh) == 0:
			print 'no good starting points within 4A of 2nd input point'
			sys.exit(2)
		gpf2 = testmesh[0] # good fill point coords
		loadpoint2 = [bestmap[gpf2[1],gpf2[2],gpf2[3]],gpf2[1],gpf2[2],gpf2[3],bestatom[gpf2[1],gpf2[2],gpf2[3]]]
		print 'second user point loaded'
		# exclude steric overlap points (use only 0 and negative affinity)
		# add points from first to second along legal pathways - then bring fill up to input volume
		spx = loadpoint1[1]
		spy = loadpoint1[2]
		spz = loadpoint1[3]
		startingpoint = (spx,spy,spz)
		flood = []
		startpoint = [bestmap[spx,spy,spz],spx,spy,spz,bestatom[spx,spy,spz]]
		flood.append(startpoint)
		# now generate neighbors - save in growing list of all neighbors
		# add neighbor point to flood that gets closer to loadpoint2
		nbor = []
		reached_B = 0
		for i in range(volume-1):
			nborfinder2(spx,nx,spy,ny,spz,nz,bestmap,bestatom,nbor,flood)
			# now sort neighbors list for min distance and add to flood list
			if reached_B == 1:
				nbor.sort()
			else:
				nborD = []
				for n in nbor:
					DtoB = sqrt((n[1] - gpf2[1])**2.0 + (n[2] - gpf2[2])**2.0 + (n[3] - gpf2[3])**2.0)
					if DtoB == 0.0:
						reached_B = 1
					nborD.append([DtoB,n[1],n[2],n[3]])
				# sort nborD and reload nbor in that order
				nborD.sort()
				nbor = []
				for n in nborD:
					node = [bestmap[n[1],n[2],n[3]],n[1],n[2],n[3],bestatom[n[1],n[2],n[3]]]
					nbor.append(node)
			# node = nbor.pop(0)
			node = nbor[0]
			Hloopcount = 0
			nborsize = len(nbor)
			while node[4] == 1:  # H atom
				x = node[1]
				y = node[2]
				z = node[3]
				Hstable = 0
				if x+1 <= nx-1:
					Hnode = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
					if Hnode in flood and bestatom[x+1,y,z] != 1:
						Hstable = 1
				if x-1 >= 0:
					Hnode = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
					if Hnode in flood and bestatom[x-1,y,z] != 1:
						Hstable = 1
				if y+1 <= ny-1:
					Hnode = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
					if Hnode in flood and bestatom[x,y+1,z] != 1:
						Hstable = 1
				if y-1 >= 0:
					Hnode = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
					if Hnode in flood and bestatom[x,y-1,z] != 1:
						Hstable = 1
				if z+1 <= nz-1:
					Hnode = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
					if Hnode in flood and bestatom[x,y,z+1] != 1:
						Hstable = 1
				if z-1 >= 0:
					Hnode = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
					if Hnode in flood and bestatom[x,y,z-1] != 1:
						Hstable = 1
				if Hstable == 0:
					Hloopcount = Hloopcount + 1
					if Hloopcount == nborsize:
						badHstartFlag = 1
						break
					node = nbor[Hloopcount]
					continue
				else:
					break
			if badHstartFlag == 1:
				break
			node = nbor.pop(Hloopcount)
			flood.append(node)
			spx = node[1]
			spy = node[2]
			spz = node[3]
		if badHstartFlag == 1:
			badHstartFlag = 0   # reset flag
			continue
		# need to finish neighbor nodes for last flood point
		nborfinder(spx,nx,spy,ny,spz,nz,bestmap,bestatom,nbor,flood)
		nbor.sort()
		# update nbor
		nbor = []
		TotalEnergy = 0.0
		for fl in flood:
			nborx = fl[1]
			nbory = fl[2]
			nborz = fl[3]
			nborfinder(nborx,nx,nbory,ny,nborz,nz,bestmap,bestatom,nbor,flood)
			TotalEnergy = TotalEnergy + fl[0]
	else:
		spx = meshpoint[1]
		spy = meshpoint[2]
		spz = meshpoint[3]
		if bestmap[spx,spy,spz] >= 0.0:
			unusablepointcount = unusablepointcount + 1
			continue
		# the start point should now be:  [spx,spy,spz]
		startingpoint = (spx,spy,spz)
		flood = []
		startpoint = [bestmap[spx,spy,spz],spx,spy,spz,bestatom[spx,spy,spz]]
		flood.append(startpoint)
		# now generate neighbors - save in growing list of all neighbors
		nbor = []
		TotalEnergy = startpoint[0]
		
		for i in range(volume-1):
			#this is where we need to build the display
			#flood is list of list x,y,z,flood[j][1:4] flood[j][4] = atom type
			# GUI
			if movie:
				save_flood(flood)
			nborfinder(spx,nx,spy,ny,spz,nz,bestmap,bestatom,nbor,flood)
			# now sort neighbors list for min energy and add to flood list
			nbor.sort()
			node = nbor[0]
			Hloopcount = 0
			nborsize = len(nbor)
			while node[4] == 1:  # H atom
				x = node[1]
				y = node[2]
				z = node[3]
				Hstable = 0
				if x+1 <= nx-1:
					Hnode = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
					if Hnode in flood and bestatom[x+1,y,z] != 1:
						Hstable = 1
				if x-1 >= 0:
					Hnode = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
					if Hnode in flood and bestatom[x-1,y,z] != 1:
						Hstable = 1
				if y+1 <= ny-1:
					Hnode = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
					if Hnode in flood and bestatom[x,y+1,z] != 1:
						Hstable = 1
				if y-1 >= 0:
					Hnode = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
					if Hnode in flood and bestatom[x,y-1,z] != 1:
						Hstable = 1
				if z+1 <= nz-1:
					Hnode = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
					if Hnode in flood and bestatom[x,y,z+1] != 1:
						Hstable = 1
				if z-1 >= 0:
					Hnode = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
					if Hnode in flood and bestatom[x,y,z-1] != 1:
						Hstable = 1
				if Hstable == 0:
					Hloopcount = Hloopcount + 1
					if Hloopcount == nborsize:
						badHstartFlag = 1
						break
					node = nbor[Hloopcount]
					continue
				else:
					break
			if badHstartFlag == 1:
				break
			node = nbor.pop(Hloopcount)
			flood.append(node)
			spx = node[1]
			spy = node[2]
			spz = node[3]
			TotalEnergy = TotalEnergy + node[0]
		if badHstartFlag == 1:
			badHstartFlag = 0   # reset flag
			continue
		# need to finish neighbor nodes for last flood point
		nborfinder(spx,nx,spy,ny,spz,nz,bestmap,bestatom,nbor,flood)
		nbor.sort()
	# now let the flooded area migrate to higher affinity area
	flagPP = 1
	PPcount = 0
	while flagPP == 1:
		# GUI
		if movie:
			save_flood(flood)
		flagPP = 0
		flood.sort()
		worstpoint = flood[volume-1]
		bestpoint = nbor[0]
		loopcount = 0
		Hloopcount = 0
		EPVloopcount = 0
		# find the energy per volume for initial fill
		TotalVolume = []
		for fl in flood:
			for CML in CoordModsList:
				point = [fl[1]+CML[0],fl[2]+CML[1],fl[3]+CML[2]]
				if point not in TotalVolume:
					TotalVolume.append(point)
		NumberVolumePoints = len(TotalVolume)
		volangstrom = NumberVolumePoints * (spacing**3)
		ThisEPV = TotalEnergy/volangstrom
		while bestpoint[0] < worstpoint[0]:
			# GUI
			if movie:
				save_flood(flood)			

			# check to see if removing worstpoint from flood breaks continuity
			#   look at state varible contig: 0 = removable, 1 = non-removable
 			contig = setstate(worstpoint,flood)
			if start_point_flag == 2:
				if worstpoint == loadpoint1:
					contig = 1
				if worstpoint == loadpoint2:
					contig = 1
			# if still contiguous: remove worstpoint and add bestpoint to flood
			# then recalculate nbor shell
			# but not if bestpoint is H and only neighbor is H atom
			if bestpoint[4] == 1 and contig == 0:
				x = bestpoint[1]
				y = bestpoint[2]
				z = bestpoint[3]
				Hstable = 0
				if x+1 <= nx-1:
					node = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
					if node in flood and bestatom[x+1,y,z] != 1:
						Hstable = 1
				if x-1 >= 0:
					node = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
					if node in flood and bestatom[x-1,y,z] != 1:
						Hstable = 1
				if y+1 <= ny-1:
					node = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
					if node in flood and bestatom[x,y+1,z] != 1:
						Hstable = 1
				if y-1 >= 0:
					node = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
					if node in flood and bestatom[x,y-1,z] != 1:
						Hstable = 1
				if z+1 <= nz-1:
					node = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
					if node in flood and bestatom[x,y,z+1] != 1:
						Hstable = 1
				if z-1 >= 0:
					node = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
					if node in flood and bestatom[x,y,z-1] != 1:
						Hstable = 1
				if Hstable == 0:
					Hloopcount = Hloopcount + 1
					bestpoint = nbor[Hloopcount]
					continue
			if contig == 0:
				# first check if bestpoint is only next to worstpoint
				x = bestpoint[1]
				y = bestpoint[2]
				z = bestpoint[3]
				Stable = 0
				if x+1 <= nx-1:
					node = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
					if node in flood:
						Stable += 1
					if node == worstpoint:
						Stable -= 1
				if x-1 >= 0:
					node = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
					if node in flood:
						Stable += 1
					if node == worstpoint:
						Stable -= 1
				if y+1 <= ny-1:
					node = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
					if node in flood:
						Stable += 1
					if node == worstpoint:
						Stable -= 1
				if y-1 >= 0:
					node = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
					if node in flood:
						Stable += 1
					if node == worstpoint:
						Stable -= 1
				if z+1 <= nz-1:
					node = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
					if node in flood:
						Stable += 1
					if node == worstpoint:
						Stable -= 1
				if z-1 >= 0:
					node = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
					if node in flood:
						Stable += 1
					if node == worstpoint:
						Stable -= 1
				if Stable <= 0:
					loopcount = loopcount + 1
					if (volume - loopcount) < 1:
						break
					worstpoint = flood[volume -1 -loopcount]
					continue
				# check if epv is beter with worstpoint changed to bestpoint
				# if not, use next best nbor point
				TotalVolume = []
				for fl in flood:
					for CML in CoordModsList:
						point = [fl[1]+CML[0],fl[2]+CML[1],fl[3]+CML[2]]
						if point not in TotalVolume:
							TotalVolume.append(point)
				NumberVolumePoints = len(TotalVolume)
				volangstrom = NumberVolumePoints * (spacing**3)
				epv = TotalEnergy/volangstrom
				if epv > ThisEPV:
					EPVloopcount = EPVloopcount + 1
					bestpoint = nbor[EPVloopcount]
					continue
				# now add bestpoint and remove worstpoint
				flood.pop(volume-1-loopcount)
				flood.append(bestpoint)
				flood.sort()
				# update nbor
				nbor = []
				for fl in flood:
					nborx = fl[1]
					nbory = fl[2]
					nborz = fl[3]
					nborfinder(nborx,nx,nbory,ny,nborz,nz,bestmap,bestatom,nbor,flood)
				nbor.sort()
				worstpoint = flood[volume-1]
				bestpoint = nbor[0]
				loopcount = 0  # reset every time a point is moved
				Hloopcount = 0 # reset every time a point is moved
			else:
				loopcount = loopcount + 1
				if (volume - loopcount) < 1:
					break
				worstpoint = flood[volume -1 -loopcount]
		# Psuedo Pod method - reach out 6A away and look for beter energy
		# only run Psuedo Pod method on fills > 11 points
		if volume <= 11:
			continue
		# only let the psuedo pods run 100 times
		if PPcount > 100:
			print "Psuedo Pod over 100 - continuing..."
			continue
		PPcount = PPcount + 1
		# first, add up energy of 10 worst points
		flood.sort()
		wplist = []
		wpvalue = []
		wp = 0
		wpcount = 0
		# create temp flood list and remove worst points so that contig will not use removed points
		tempwpflood = []
		for fl in flood:
			tempwpflood.append(fl)
		while wpcount < volume:
			worstpoint = flood[volume -1 -wpcount]
			contig = setstate(worstpoint,tempwpflood)
			if start_point_flag == 2:
				if worstpoint == loadpoint1:
					contig = 1
				if worstpoint == loadpoint2:
					contig = 1
			if contig == 0:
				tempwpflood.remove(worstpoint)
				if wp == 0:
					wp1 = worstpoint[0]
					wpvalue.append(wp1)
					wplist.append(worstpoint)
				if wp == 1:
					wp2 = worstpoint[0] + wp1
					wpvalue.append(wp2)
					wplist.append(worstpoint)
				if wp == 2:
					wp3 = worstpoint[0] + wp2
					wpvalue.append(wp3)
					wplist.append(worstpoint)
				if wp == 3:
					wp4 = worstpoint[0] + wp3
					wpvalue.append(wp4)
					wplist.append(worstpoint)
				if wp == 4:
					wp5 = worstpoint[0] + wp4
					wpvalue.append(wp5)
					wplist.append(worstpoint)
				if wp == 5:
					wp6 = worstpoint[0] + wp5
					wpvalue.append(wp6)
					wplist.append(worstpoint)
				if wp == 6:
					wp7 = worstpoint[0] + wp6
					wpvalue.append(wp7)
					wplist.append(worstpoint)
				if wp == 7:
					wp8 = worstpoint[0] + wp7
					wpvalue.append(wp8)
					wplist.append(worstpoint)
				if wp == 8:
					wp9 = worstpoint[0] + wp8
					wpvalue.append(wp9)
					wplist.append(worstpoint)
				if wp == 9:
					wp10 = worstpoint[0] + wp9
					wpvalue.append(wp10)
					wplist.append(worstpoint)
					wp = 10
					break
				wp = wp + 1
			wpcount = wpcount + 1
		# next, loop through flood points and project pseudopods out to 6A
		# but not from the worst 10 nodes or from H atoms
		bestPP = []
		for fl in flood:
			if fl in wplist:
				continue
			if fl[4] == 1:   # atomtype 1 is an H atom
				continue
			totenergyPP = 0.0
			lengthPP = 0
			PPlist = []
			# generate 6A shell around point, exclude points in protein and flood
			Shell = []
			for CML6 in CoordModsList6A:
				x = fl[1]+CML6[0]
				if x >= nx or x < 0:
					continue
				y = fl[2]+CML6[1]
				if y >= ny or y < 0:
					continue
				z = fl[3]+CML6[2]
				if z >= nz or z < 0:
					continue
				node = [bestmap[x,y,z],x,y,z,bestatom[x,y,z]]
				if node[0] > 0.0:
					continue
				if node in flood:
					continue
				Shell.append(node)
			# pick best point in shell
			Shell.sort()
			bestShellpoint = Shell[0]
			# march from flood point to best Shell point - record length and total energy
			currentPPnode = [fl[0],fl[1],fl[2],fl[3],fl[4]]
			bspx = bestShellpoint[1]
			bspy = bestShellpoint[2]
			bspz = bestShellpoint[3]
			for i in range(10):
				x = currentPPnode[1]
				y = currentPPnode[2]
				z = currentPPnode[3]
				# check each direction from current point to get closer to best point
				checkNode = []
				subchecknode = []
				if x+1 <= nx-1:
					nodePP = [bestmap[x+1,y,z],x+1,y,z,bestatom[x+1,y,z]]
					if nodePP in Shell:
						DtoBSPsq = (nodePP[1] - bspx)**2.0 + (nodePP[2] -bspy)**2.0 + (nodePP[3] - bspz)**2.0
						subchecknode = [DtoBSPsq,nodePP]
						checkNode.append(subchecknode)
				if x-1 >= 0:
					nodePP = [bestmap[x-1,y,z],x-1,y,z,bestatom[x-1,y,z]]
					if nodePP in Shell:
						DtoBSPsq = (nodePP[1] - bspx)**2.0 + (nodePP[2] -bspy)**2.0 + (nodePP[3] - bspz)**2.0
						subchecknode = [DtoBSPsq,nodePP]
						checkNode.append(subchecknode)
				if y+1 <= ny-1:
					nodePP = [bestmap[x,y+1,z],x,y+1,z,bestatom[x,y+1,z]]
					if nodePP in Shell:
						DtoBSPsq = (nodePP[1] - bspx)**2.0 + (nodePP[2] -bspy)**2.0 + (nodePP[3] - bspz)**2.0
						subchecknode = [DtoBSPsq,nodePP]
						checkNode.append(subchecknode)
				if y-1 >= 0:
					nodePP = [bestmap[x,y-1,z],x,y-1,z,bestatom[x,y-1,z]]
					if nodePP in Shell:
						DtoBSPsq = (nodePP[1] - bspx)**2.0 + (nodePP[2] -bspy)**2.0 + (nodePP[3] - bspz)**2.0
						subchecknode = [DtoBSPsq,nodePP]
						checkNode.append(subchecknode)
				if z+1 <= nz-1:
					nodePP = [bestmap[x,y,z+1],x,y,z+1,bestatom[x,y,z+1]]
					if nodePP in Shell:
						DtoBSPsq = (nodePP[1] - bspx)**2.0 + (nodePP[2] -bspy)**2.0 + (nodePP[3] - bspz)**2.0
						subchecknode = [DtoBSPsq,nodePP]
						checkNode.append(subchecknode)
				if z-1 >= 0:
					nodePP = [bestmap[x,y,z-1],x,y,z-1,bestatom[x,y,z-1]]
					if nodePP in Shell:
						DtoBSPsq = (nodePP[1] - bspx)**2.0 + (nodePP[2] -bspy)**2.0 + (nodePP[3] - bspz)**2.0
						subchecknode = [DtoBSPsq,nodePP]
						checkNode.append(subchecknode)
				if len(checkNode) == 0:
					break
				lengthPP = lengthPP + 1
				checkNode.sort()
				subchecknode = checkNode[0]
				nodePP = subchecknode[1]
				PPlist.append(nodePP)
				# remove nodePP from Shell
				Shell.remove(nodePP)
				totenergyPP = totenergyPP + nodePP[0]
				if subchecknode[0] <= 0.01:
					# reached the best node
					break
				# stop at the first H reached
				if nodePP[4] == 1:
					break
				currentPPnode = [nodePP[0],nodePP[1],nodePP[2],nodePP[3],nodePP[4]]
			bestnodePP = [totenergyPP,PPlist,lengthPP]
			bestPP.append(bestnodePP)
		bestPP.sort()
		if len(bestPP) == 0:
			break
		bestnodePP = bestPP.pop(0)
		# if the pseudopod is better than up to 10 worst, remove worst and add pseudopod
		if bestnodePP[2] > len(wpvalue):
			break
		if bestnodePP[0] < wpvalue[bestnodePP[2]-1]:
			for i in range(bestnodePP[2]):
				flood.remove(wplist[i])
			for nodePP in bestnodePP[1]:
				flood.append(nodePP)
			# now update nbor and start over to migrate
			nbor = []
			for fl in flood:
				nborx = fl[1]
				nbory = fl[2]
				nborz = fl[3]
				nborfinder(nborx,nx,nbory,ny,nborz,nz,bestmap,bestatom,nbor,flood)
			nbor.sort()
			flagPP = 1 # repeat the migration loop to update with the pseudopod

	# now output flood as xyz file
	i = 0
	TotalEnergy = 0.0
	for fl in flood:
		TotalEnergy = TotalEnergy + fl[0]
	spconvx = (startingpoint[0] - xcent)*spacing + centerx
	spconvy = (startingpoint[1] - ycent)*spacing + centery
	spconvz = (startingpoint[2] - zcent)*spacing + centerz
	PDBsublist = []
	for fl in flood:
		i = i + 1
		x = (fl[1] - xcent)*spacing + centerx
		y = (fl[2] - ycent)*spacing + centery
		z = (fl[3] - zcent)*spacing + centerz
		if fl[4] == 7:
			atomchr = 'N'
		if fl[4] == 6:
			atomchr = 'S'
		if fl[4] == 5:
			atomchr = 'A'
		if fl[4] == 4:
			atomchr = 'O'
		if fl[4] == 3:
			atomchr = 'P'
			# note, this will color the N atom pink (the PDB color for Phosphorus)
		if fl[4] == 2:
			atomchr = 'C'
		if fl[4] == 1:
			atomchr = 'H'
		afin = abs(fl[0])
		PDBsublist.append([i,atomchr,x,y,z,afin])
	PDBlist.append(PDBsublist)
	# now calculate volume of flood with 1A radii
	TotalVolume = []
	for fl in flood:
		for CML in CoordModsList:
			point = [fl[1]+CML[0],fl[2]+CML[1],fl[3]+CML[2]]
			if point not in TotalVolume:
				TotalVolume.append(point)
	NumberVolumePoints = len(TotalVolume)
	volangstrom = NumberVolumePoints * (spacing**3)
	epv = TotalEnergy/volangstrom
	OUTcount = OUTcount + 1
	OUTlist.append([epv,volangstrom,OUTcount])

# GUI
if movie:
	save_flood(flood)
	output.close()

# sort and generate outputs
OUTlist.sort()
OUTlength = len(OUTlist)
lastOP = OUTlist[0]
# find non-duplicates
newOUTlist = []
newOUTlist.append(lastOP)
for i in range(OUTlength - 1):
	x = i + 1
	thisOP = OUTlist[x]
	if thisOP[0] == lastOP[0]:
		continue
	else:
		lastOP = thisOP
		newOUTlist.append(thisOP)
resultName = FileBaseName + '_' + str(volume) + 'Results.txt'
results_out = open(resultName,'w')
i = 0
if start_point_flag == 0:
	print '10 Unique Results:  (note, some starting points converge to same answer, so may be less than 10)'
for OP in newOUTlist:
	print 'Output #',i+1,' Total Volume = ',OP[1],' Total Energy per Vol, EPV = ',OP[0],' (Kcal/mol A**3)'
	results_out.write('Output #%3i, Total Volume = %7.2f, Total Energy per Vol, EPV = %15.8f\n' 
			   % (i+1,OP[1],OP[0]))
	# point loops
	PDBname = 'FILL_' + str(volume) + 'out' + str(i+1) + '.pdb'
	PDB_out = open(PDBname,'w')
	for PDB in PDBlist[OP[2]-1]:
		PDB_out.write('ATOM%7d  %s   FIL     1%12.3f%8.3f%8.3f  1.0%8.3f     0.001 %s\n' %
				(PDB[0],PDB[1],PDB[2],PDB[3],PDB[4],PDB[5],PDB[1]))
	PDB_out.close()
	i = i + 1
# clean house
results_out.close()
if Sflag == 1: Smap.close()
if Aflag == 1: Amap.close()
if Cflag == 1: Cmap.close()
if Hflag == 1: Hmap.close()
if Nflag == 1: Nmap.close()
if NAflag == 1: NAmap.close()
if Oflag == 1: Omap.close()
emap.close()
dmap.close()
