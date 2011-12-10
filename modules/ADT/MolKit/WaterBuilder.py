############################################################################
#
# Authors: Stefano Forli, Ruth Huey
#
# Copyright: M. Sanner TSRI 2010
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/WaterBuilder.py,v 1.1 2010/10/07 17:28:48 rhuey Exp $
#
# $Id: WaterBuilder.py,v 1.1 2010/10/07 17:28:48 rhuey Exp $
#

"""
This module implements the WaterBuilder class which adds water atom records to pdb files 

"""

from numpy import sin, cos, sqrt, array, sum, cross, radians, dot
import math, os


class WaterBuilder:
    """Class for adding waters to a pdb file.
    """



    def __init__(self, space=3.0, bond_dist=1.85, pcharge=0.000, residue="WAT", ATYPE="W", GPF=False, verbose=False):
        self.space = space
        self.bond_dist = bond_dist
        self.pcharge = pcharge
        self.residue = residue
        self.ATYPE = ATYPE
        self.GPF = GPF
        self.verbose = verbose


    def atom_coord(self, atom):
        return map(float, atom[28:56].split())


    def dist(self, firstline, secondline, precision=4):  
        coord1, coord2 = map(self.atom_coord, [firstline, secondline])
        return round(sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2), precision)


    def mean_pdb(self, firstline, secondline):
        # INFO   : calculate the mean point between two PDB atom lines
        # INPUT  : two pdb lines, an number for residue and atom numbering
        # OUTPUT : a pdb line
        coord1 = self.atom_coord(firstline)
        coord2 = self.atom_coord(secondline)
        x = (coord1[0]+coord2[0])/2
        y = (coord1[1]+coord2[1])/2
        z = (coord1[2]+coord2[2])/2
        atype = firstline[12:16]
        residue = "MEA"
        self.chain = "Y"
        count = 1
        index = 1
        #?SHOULD THIS BE atype from line 58?
        mean_atom = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %1s" % (count, atype, residue, self.chain, index, x, y, z, self.ATYPE)
        #mean_atom = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %1s" % (count, self.ATYPE, residue, self.chain, index, x, y, z, self.ATYPE)
        return mean_atom


    def closest(self, first_atom, atom_list, cutoff=99999):
        # INFO   : find the self.closest atom 
        # INPUT  : a PDB atom, a list of PDB atoms [, cutoff self.distance]
        # OUTPUT : the self.closest atom [=null, if cutoff not satisfied] and short self.distance found
        # EXTRA  : self.dist function required
        best_self.distance = 999999
        best_candidate = None
        for second_atom in atom_list:
            self.distance = self.dist(first_atom, second_atom)
            if self.distance < best_self.distance:
                best_self.distance = self.distance
                if best_self.distance < cutoff:
                    best_candidate = second_atom
        #if best_candidate != "":
        #    return best_candidate, best_self.distance
        #else:
        #    return best_candidate, best_self.distance
        return best_candidate, best_self.distance


    def rotatePoint(self, pt,m,ax):
        # From Ludo
        x = pt[0]
        y = pt[1]
        z = pt[2]
        u = ax[0]
        v = ax[1]
        w = ax[2]
        ux = u*x
        uy = u*y
        uz = u*z
        vx = v*x
        vy = v*y
        vz = v*z
        wx = w*x
        wy = w*y
        wz = w*z
        sa = sin(ax[3])
        ca = cos(ax[3])
        pt[0] = (u*(ux+vy+wz)+(x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa) +m[0]
        pt[1] = (v*(ux+vy+wz)+(y*(u*u+w*w)-v*(ux+wz))*ca+(wx-uz)*sa)  +m[1]
        pt[2] = (w*(ux+vy+wz)+(z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa) +m[2]
        return pt


    def mean3(self, firstline, secondline, thirdline):
        # INFO   : calculate the mean point between two PDB atom lines
        # INPUT  : two pdb lines, an number for residue and atom numbering
        # OUTPUT : a pdb line
        coord1 = self.atom_coord(firstline)
        coord2 = self.atom_coord(secondline)
        coord3 = self.atom_coord(thirdline)
        x = (coord1[0]+coord2[0]+coord3[0])/3
        y = (coord1[1]+coord2[1]+coord3[1])/3
        z = (coord1[2]+coord2[2]+coord3[2])/3
        atype = firstline[12:16]
        residue = "MEA"
        chain = "Y"
        count = 1
        index = 1
        mean_atom = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %1s" % (count, self.ATYPE, residue, chain, index, x, y, z, self.ATYPE)
        return mean_atom


    def hydro(self, atom1, atom2, atom3=None, spacing=3):
        # INFO   : place a water atom W at given self.distance from an atom
        # INPUT  : (a) two pdb lines for N-H; (b) three pdb lines for -[N|O]- acceptor, -O-H donors; the spacing self.distance between the water and the atom
        # OUTPUT : one pdb line 
        #
        # Synopsis:
        #
        # -C-NH-C  => atom1 = N, atom2 = H
        # -C-N=C   => atom1 = C, atom2 = N, atom3 = C
        # -C-O-H   => atom1 = C, atom2 = O, atom3 = H

        #if not atom3:
        #    atype = line.split()[-1]
        if atom2.split()[-1]=="HD":
            spacing -= 1 # to put the W at 3A from the N-H
        index = 99 
        coord2 = self.atom_coord(atom2)
        x =(coord2[0])
        y =(coord2[1])
        z =(coord2[2])
        if self.verbose: print x, y, z
        #atype = "W"
        #residue="WAT"
        # residue 17:20
        # chain 21
        chain = atom1[21]
        residue = atom1[17:20]
        #chain="X"
        if atom3:
            atom4 = self.mean_pdb(atom1, atom3)
            vec_module = self.dist(atom2, atom4)
            coord1 = self.atom_coord(atom4) 
        else:
            coord1 = self.atom_coord(atom1)
            vec_module = self.dist(atom1, atom2)
        alpha = math.acos((coord2[0]-coord1[0])/vec_module)    # x-axis angle
        beta  = math.acos((coord2[1]-coord1[1])/vec_module)    # y-axis angle
        gamma = math.acos((coord2[2]-coord1[2])/vec_module)    # z-axis angle
        wat_x = spacing*cos(alpha)+x
        wat_y = spacing*cos(beta)+y
        wat_z = spacing*cos(gamma)+z
        wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, index, self.ATYPE, residue, chain, index, wat_x, wat_y, wat_z, self.pcharge, self.ATYPE)
        #wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, index, self.ATYPE, self.residue, chain, index, wat_x, wat_y, wat_z, self.pcharge, self.ATYPE)
        return wet


    def hydroH(self, atom1, atom2, atom3):
        middle = self.hydro(atom1, atom2, spacing = 0)
        avg = self.mean_pdb(middle, atom3)
        last = self.hydro(avg, atom2, spacing=self.space)
        if self.verbose:
            print "middle=", middle
            print "avg=", avg
            print "last=", last


    def vector(self, p1, p2=None):
        # accept both atoms and self.vectors as input
        #
        #if p1 == 0:
        #    coord1 = [ 0.0, 0.0, 0.0 ]
        #else:
        if self.verbose: print "received ", p1, p2
        if type(p1) == type(str()):
            p1 = self.atom_coord(p1)
        x1,y1,z1 = p1
        if type(p2) == type(str()):
            p2 = self.atom_coord(p2)
        if not p2 == None:
            x2,y2,z2 = p2
            vec_x = x2-x1
            vec_y = y2-y1
            vec_z = z2-z1
            # it must be an array
            vec = array([vec_x, vec_y, vec_z], 'f')
            if self.verbose: print "REAL VECTOR", vec
        else:
            vec = array([p1[0], p1[1], p1[2] ], 'f' )
            if self.verbose: print "ATOM VECTOR", vec
        return vec


    def norm(self, A):
        "Return vector norm"
        return sqrt(sum(A*A))


    def normalize(self, A):
        "Normalize the Vector"
        return A/self.norm(A)


    def bound(self, atom, structure, exclude=None):
        """
        Identify all the atoms in "structure" that are @bond_dist from "atom"
        """
        bound_list = []
        tolerance = 0
        bond_dist = self.bond_dist
        if atom.split()[-1] == "HD":
            bond_dist = 1.15 # previous bond self.dist of 1.1 wouldn't be big enough for -S-H
            if self.verbose: print "HD mode"
        for candidate in structure:
            if candidate==atom or candidate==exclude:
                pass
            else:
                if candidate[0:4]=="ATOM" or candidate[0:6]=="HETATM":
                    if candidate.split()[-1]=="SA" or candidate.split()[-1]=="S":
                        tolerance = .35
                    else:
                        tolerance = 0
                    if self.dist(atom, candidate)<=bond_dist+tolerance:
                        if not candidate in bound_list:
                            bound_list.append(candidate)
                    else:
                        pass
        if len(bound_list) > 0:
            return bound_list
        else:
            print "ERROR: this atom seems to be disconnected:"
            print atom
            print "exit 3"
            exit(1)


    def calc_plane(self, atom1, atom2, atom3):
        # weird but it works...
        v12 = self.vector(atom1, atom2)
        v13 = self.vector(atom3, atom2)
        plane = cross(v12, v13)
        plane = self.normalize(plane)
        residue = atom1[17:20]
        chain = atom1[21]
        #plane = plane/(sqrt(sum(plane*plane)))
        if self.verbose:
            print atom1
            print atom2
            print atom3
            print "PLANE FREE> coords: ", plane
            print "PLANE FREE> type: ", type(plane)
            print "PLANE FREE> atoms:"
            self.keyw, index, self.ATYPE, residue, chain, self.pcharge = "ATOM  ", 1, "C", residue, chain, 0.0000
            #self.keyw, index, self.ATYPE, self.residue, self.chain, self.pcharge = "ATOM  ", 1, "C", "RES", 1, 0.0000
            print "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s" % (self.keyw, index, self.ATYPE, residue, chain, index, plane[0], plane[1], plane[2], self.pcharge, self.ATYPE)
            #print "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s" % (self.keyw, index, self.ATYPE, self.residue, self.chain, index, plane[0], plane[1], plane[2], self.pcharge, self.ATYPE)
            print "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, index, "X" , residue, chain, index, 0,0,0 , self.pcharge, "X")

            centroid = self.mean3(atom1, atom2, atom3)
            arrow = self.vector(atom1, centroid)+ plane
            #arrow = vec_sum(self.vector(atom1, centroid), plane)
            print "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % ("ATOM  ", 1, "A", 1, 1, 1, arrow[0], arrow[1], arrow[2], 0.000, "A")
        return plane


    def vec_sum(self, vec1, vec2):
        return array([vec1[0]+vec2[0], vec1[1]+vec2[1], vec1[2]+vec2[2] ], 'f')


    def coplanar(self, plane, structure, reference, tolerance=.2):
        # identify all the atoms in "structure" that are co-planar with "plane"
        # the dot product between the point and the plane must be ~= 0
        coplane_list = []
        for atom in structure:
            position = self.vector(reference, atom)
            if self.verbose:
                print "POSITION", position
                print self.atom_coord(atom)
            if dot(plane, position) <= tolerance:
                coplane_list.append(atom)
        return coplane_list


    #def dot(self.vector1, self.vector2):
    #    dot_product = 0.
    #    for i in range(0, len(self.vector1)):
    #        dot_product += (self.vector1[i] * self.vector2[i])
    #    return dot_product


    def furanbolic(self, atom, structure, max=2.35):
        # it should walk with pre-filtered co-planar atoms
        # HD's are automatically excluded
        the_ring = [atom]
        #print "I WAS REQUESTED TO OPERATE ON ", atom
        #print "I WILL SEARCH ON %d ATOMS"% len(structure)
        for item in structure:
            if item!=atom:
                if item.split()[-1]!="HD":
                    if self.dist(atom, item)<max:
                        the_ring.append(item)
                        if self.verbose:
                            print "APPEND"
                    #else:
                        #print "failed self.distance"
                        #print item
        #print "\n ########## FURANBOLIC CHECKPOINT ### \n"
        #print the_ring
        #print len(the_ring)
        if len(the_ring)==5:
            print " - possible furan/oxazole found..."
            return True
        if len(the_ring)>6:
            print "WARNING: multiple atoms match the furan/oxazole check..."
            return True
        else:
            return False


    def Osp2(self, oxygen, atom1, atom2):
        waters = []
        # self.hydroxyl/ether mode
        angles=[120, -120]
        #angles = range(0, 360, 10)
        oxyvector = self.vector(oxygen, atom1)
        oxyvector = self.normalize(oxyvector)
        residue = atom1[17:20]
        chain = atom1[21]
        for a in angles:
            roto = [oxyvector[0], oxyvector[1], oxyvector[2],  radians(a)]
            lone_pair_vector = self.vector(atom2, oxygen)
            lone_pair_vector = self.normalize(lone_pair_vector)
            water =  self.rotatePoint(-lone_pair_vector*self.space, self.atom_coord(oxygen), roto)  
            #residue = "99"
            #chain = "1"
            wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n"%(self.keyw, 1, self.ATYPE, residue, chain, 1, water[0], water[1], water[2], self.pcharge, self.ATYPE)
            #wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n"%(self.keyw, 1, self.ATYPE, residue, chain, 1, water[0], water[1], water[2], self.pcharge, self.ATYPE)
            waters.append(wet)
        if self.verbose: print "osp2 returning %d waters"%(len(waters)) 
        return waters


    def Osp2_NEW(self, oxygen, atom1, atom2):
        waters = []
        # self.hydroxyl/ether mode
        #angles = [120, -120]
        angles = range(0, 360, 10)
        oxyvector = self.vector(oxygen, atom1)
        oxyvector = self.vector(atom1, atom2)
        oxyvector = self.normalize(oxyvector)
        mid = self.mean_pdb(atom1, atom2)
        residue = atom1[17:20]
        chain = atom1[21]
        for a in angles:
            roto = [oxyvector[0], oxyvector[1], oxyvector[2],  radians(a)]
            lone_pair_vector = self.vector(mid, oxygen)
            lone_pair_vector = self.normalize(lone_pair_vector)
            water =  self.rotatePoint(+lone_pair_vector*self.space, self.atom_coord(oxygen), roto)  
            #residue = "99"
            #chain = "1"
            # change to self.residue? replacing 99 by WAT??
            wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n"%(self.keyw, 1, self.ATYPE, residue, chain, 1, water[0], water[1], water[2], self.pcharge, self.ATYPE)
            waters.append(wet)
        if self.verbose: print "Osp2_NEW returning %d waters"%(len(waters)) 
        return waters


    def gpfminmax(self, gpf):
        # INPUT : gpf_file
        # OUTPUT: box coordinates (x,y,z), (X,Y,Z)
        """ Return max/min coordinates of
        the box described in the GPF"""
        if not self.GPF:
            return True
        if self.verbose: print 'gpf=', gpf    
        file = open(gpf, 'r')
        lines = file.readlines()
        file.close()
        for line in lines:
            tmp=line.split()
            #print tmp
            if tmp[0] == "gridcenter":
                center_x, center_y, center_z = map(float, tmp[1:4])
            if tmp[0] == "npts":
                pts_x, pts_y, pts_z = map(int, tmp[1:4])
            if tmp[0] == "spacing":
                res = float(tmp[1])
                if self.verbose: print "spacing: ", res
        step_x = pts_x/2 * res
        step_y = pts_y/2 * res
        step_z = pts_z/2 * res
        x_min = center_x - step_x
        x_max = center_x + step_x
        y_min = center_y - step_y
        y_max = center_y + step_y
        z_min = center_z - step_z
        z_max = center_z + step_z
        print " - using the GPF box filter [ %s ]"% gpf
        return [x_min, y_min, z_min], [x_max, y_max, z_max]


    def in_the_box(self, atom_list, MIN, MAX):
        # INPUT : pdb-like atom, MIN = [x,y,z], MAX = [x,y,z]
        # OUTPUT: True/False
        good = []
        for atom in atom_list:
            pos = self.atom_coord(atom)
            if MIN[0] < pos[0] < MAX[0]:
                if MIN[1] < pos[1] < MAX[1]:
                    if MIN[2] < pos[2] < MAX[2]:
                        good.append(atom)
        return good


    def hydrate(self, pdbqt, gpf=None, force=False, extended_atoms=False, output=None):
        added_waters = []
        atoms_list = []
        hydrate_list = []
        numbering_stuff = []
        water_mates = []
        self.FORCE = force
        self.EXTENDED_ATOMS = extended_atoms
        if not os.path.exists(pdbqt):
            print "\n\n\t# ERROR #\n\t Input file not found!"
        #Initialize name here:
        name = os.path.splitext(os.path.basename(pdbqt))[0]
        if gpf is not None:
            self.MIN, self.MAX = self.gpfminmax(gpf)
        else:
            if self.GPF:
                print "gpf file must be specified"
                return 
        if output is not None:
            if os.path.isdir(output):
                if self.verbose: print " - saving the file in the path => %s" % output
                name = output+os.path.sep+name+"_HYDRO.pdbqt"
            else:
                if self.verbose: print " - saving the file => %s" % output
                name = output

        input = open(pdbqt, 'r').readlines()
        for line in input:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                atype = line.split()[-1]
                atoms_list.append(line)
                if atype == "OA" or atype == "NA" or atype == "HD":
                    hydrate_list.append(line)
            if line[0:4] == "ATOM":
                self.keyw = "ATOM  "
            if line[0:6] == "HETATM":
                self.keyw = "HETATM"
        
        #9/24/2010 @@ START HERE=>
            #print " - total number of atoms : ", len(atoms_list), 
        if len(hydrate_list):
            if self.GPF:
                hydrate_list = self.in_the_box(hydrate_list, MIN, MAX)
            if self.verbose:
                print "=========="
                print " - hydratable atoms : %d / %d " % (len(hydrate_list), len(atoms_list))
                print "=========="
                for i in hydrate_list:
                    print i[:-1]
                print "=========="
        else:
            if not FORCE:
                print " [ No atoms to hydrate ]"
                #exit(0)
            else:
                print " [ No atoms to hydrate... FORCING TO SAVE OUTPUT...  ]"
            return

        # Scan the list to add waters
        for atom in hydrate_list:
            atype = atom.split()[-1]
            if self.verbose:
                print "PROCESSING ATOM :", atom
                print "atype = ", atype
            HYDROXYL = False
            # ordinal position in the original file
            position = int(atom.split()[1])
            # add the atom to be hydrated to the buffer list
            waters_generated = [atom]  #@@reinitialize waters_generated here
            if self.verbose:
                print "processing ", atom
                print "0: water_mates=", water_mates
                print "1: len(waters_generated)=", len(waters_generated)
                print "2:  len(water_mates=", len(water_mates)
            # find the master atom(s)
            master = self.bound(atom, atoms_list)
            if len(master) == 0:
                print "\n\nERROR: this atom is disconnected:\n", atom
                #??return() instead??@@
                print "exit 0"
                exit(1)
            # HYDROGENS #####################################################################################
            #print "before HD , len(master)=", len(master), " atype=", atype
            if atype == "HD":
                # check for errors
                if len(master) > 1:
                    print "\n\nERROR (HD) : there is a proximity error and the following hydrogen is in close contact with more than one atom"
                    print atom
                    print "Bound mates:"
                    for m in master:
                        print m[:-1]," ==>", self.dist(m,atom)
                    #??return() instead??@@
                    print "exit 1"
                    exit(1)
                else:
                    # calculate the Water vector
                    wet = self.hydro(master[0], atom) #@@REPAIR THIS apparent duplicate!!!
                    if self.verbose: 
                        print "after call to hydro, wet=",wet
                    waters_generated.append(wet)
                    if self.verbose:
                        print "HD: added wet: now %d waters_generated:" %len(waters_generated)
                        for ind in waters_generated: print ind
                        print "HD: added w_g to w_m: previously %d water_mates:" %len(water_mates)
                        for ind in water_mates: print ind
                        print "HD: 2: %d waters_generated:" %len(waters_generated)
                        for ind in waters_generated: print ind
                water_mates.append(waters_generated)
                if self.verbose: print "HD: 2: NOW %d water_mates:" %len(water_mates)
                #numbering_stuff has an entry for for each entry in water_mates
                numbering_stuff.append([position, (len(waters_generated)-1)] )

            # OXYGENS #####################################################################################
            if atype == "OA":
                # OA includes the following options:
                #
                ## Two mates
                #  ---------
                #
                # -X-OA-X-    (ethers)
                #
                # -X-OA-HD    (hydroxyls)
                #
                # [-X-OA-X-]  (furan/pyran like)
                #
                #
                ## One mate
                # ----------
                #
                # -X=O        (carbonyl, carboxylate, nitro)   [ DONE]
                #
                #
                #  |
                # -X=O        (sulphonyl, phosphonyl) 
                #  |
                #
                if len(master)==1:
                    # # identify the mates of the master
                    MASTER = master[0]
                    residue = MASTER[17:20]
                    chain = MASTER[21]
                    mates = self.bound(MASTER, atoms_list, exclude=atom)
                    #print "masters", master
                    #print "mates", mates
                    # Phosphates check
                    if len(mates) <=2:
                        # 1. calculate the plane
                        v12 = self.vector(atom, MASTER)
                        v23 = self.vector(mates[0],MASTER) 
                        #plane = cross(v12, v23)
                        plane0,plane1,plane2 = self.normalize(cross(v12,v23))
                        self.chain = 1 # TODO read this from the atom
                        #keyw, index, ATYPE, residue, chain, pcharge = "ATOM  ", 1, "C", "RES", 1, 0.0000
                        #print "CARBONYL"
                        #print "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s" % (keyw, index, ATYPE, residue, chain, index, plane[0], plane[1], plane[2], pcharge, ATYPE)

                        # O-lone pair 1
                        roto = [ plane0, plane1, plane2, radians(50) ]
                        wat = self.rotatePoint(self.normalize(-v12)*self.space, self.atom_coord(atom), roto)
                        wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, 1, self.ATYPE, residue, chain, 1, wat[0], wat[1], wat[2], self.pcharge, self.ATYPE)
                        #wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, 1, self.ATYPE, self.residue, self.chain, 1, wat[0], wat[1], wat[2], self.pcharge, self.ATYPE)
                        if self.verbose:
                            print "O1: wet:" 
                            print wet
                        waters_generated.append(wet)
                        if self.verbose:
                            print "waters_generated:"
                            print waters_generated
                            print "O1: waters_generated=", 
                            for ind in waters_generated: print ind
            
                        # O-lone pair 2
                        roto = [ plane0, plane1, plane2, radians(-50) ]
                        if self.verbose: print "calling rotatePoint with self.atom_coord(", atom, "), roto=",roto
                        wat = self.rotatePoint(self.normalize(-v12)*self.space, self.atom_coord(atom), roto)
                        if self.verbose: print "wat =", wat
                        wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, 1, self.ATYPE, residue, chain, 1, wat[0], wat[1], wat[2], self.pcharge, self.ATYPE)
                        #wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, 1, self.ATYPE, self.residue, self.chain, 1, wat[0], wat[1], wat[2], self.pcharge, self.ATYPE)
                        waters_generated.append(wet)
                        if self.verbose:
                            print "O2: %d wet: %s" %(len(wet), wet)
                            print "O2: %d waters_generated="%(len(waters_generated)), 
                            for ind in waters_generated: print ind
                        water_mates.append(waters_generated)
                        if self.verbose:
                            print "O2: %d water_mates="%(len(water_mates)), 
                            for ind in water_mates: print ind
                        numbering_stuff.append([int(position), len(waters_generated)-1])
                        #print waters_generated
                    else:
                        if self.EXTENDED_ATOMS:
                            #if self.verbose:
                                #print "master", master[0]
                                #print "mates", mates
                                #print atom
                                #print "Possible phosphate/sulphate atom... I can't stand Phosphates... GIVING UP"
                                #print "\n\n[ we're not equipped to manage those atoms... yet ]"
                                #print "\n\n[ we're testing to manage those atoms...EXPERIMENTAL!  ]"
                            #chain = 1
                            #self.residue = 1 
                            #print "\n\n  ############## TESTING FACILITY #########"
                            directive = self.vector(master[0],atom)
                            directive = self.normalize(directive)*self.space
                            #if self.verbose: print "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (keyw, 1, self.ATYPE, residue, chain, 1, directive[0], directive[1], directive[2], pcharge, self.ATYPE)

                            for q in mates:
                                #print q
                                position_q = self.vector(q) 
                                position_q = self.vec_sum(position_q, directive)
                                push = self.normalize(self.vector(atom, position_q))
                                start = self.vector(atom)
                                lpair = self.vec_sum(start, push*self.space)
                                wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, 1, self.ATYPE, residue, chain, 1, lpair[0], lpair[1], lpair[2], self.pcharge, self.ATYPE)
                                #wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (self.keyw, 1, self.ATYPE, self.residue, self.chain, 1, lpair[0], lpair[1], lpair[2], self.pcharge, self.ATYPE)
                                if self.verbose:
                                    print "EA: wet=", 
                                    for ind in wet: print ind
                                #position_q = normalize(position_q)
                                #ACTUAL = normalize(-vector(atom, position_q))*space
                                waters_generated.append(wet)
                                if self.verbose:
                                    print "EA: waters_generated=", 
                                    for ind in waters_generated: print ind
                                #if self.verbose: print "EXT ATOMS: waters_generated=", waters_generated
                            water_mates.append(waters_generated)
                            if self.verbose: 
                                print "661:EA: %d water_mates:"%(len(water_mates)), 
                                for ind in water_mates: 
                                    for j in ind:
                                        print j
                            #if self.verbose: print "         : water_mates=", water_mates
                            numbering_stuff.append([int(position), len(waters_generated)-1])
                            #print "\n\n  ############## TESTING FACILITY #########"
                            #exit(1)
                if len(master) == 2:
                    #print "---[ MODE 2 ]---"
                    for m in master:
                        if m.split()[-1] == "HD":
                            HYDROXYL = True
                            #print " >>> Found hydroxyl"
                    if not HYDROXYL:
                        # calculate the plane formed by oxygen and the two masters
                        O_plane = self.calc_plane(atom, master[0], master[1])
                        # get the list of all the co-planar atoms in the structure
                        coplanar_mates = self.coplanar(O_plane, atoms_list, atom)
                        #print "COPLANAR_MATES", len(coplanar_mates)
                        #print "FURANBOLIC", furanbolic(atom, coplanar_mates)
                        # check if there are at least 4 coplanar atoms to make a ring with the OA
                        if len(coplanar_mates) >=4 and self.furanbolic(atom, coplanar_mates):    
                            wet = self.hydro(master[0], atom, master[1]) #2: @@REPAIR THIS apparent duplicate!!!
                            if self.verbose:
                                print "notHYDROXYL: wet=", 
                                for ind in wet: print ind
                            waters_generated.append(wet)
                            if self.verbose:
                                print "notHYDROXYL: waters_generated=", 
                                for ind in waters_generated: print ind
                            #if self.verbose:
                            #    print "NOT_HYDROXYL: waters_generated=", waters_generated
                            #print "FURAN MODE"
                        else:
                            if self.verbose: print "FURAN MODE FAILED"
                            lp_waters = self.Osp2(atom, master[0], master[1])
                            if self.verbose:
                                print "lp_waters:",
                                for ind in lp_waters: print ind
                            for w in lp_waters:
                                waters_generated.append(w)
                            if self.verbose:
                                print "else in not HYDROXYL: waters_generated=", 
                                for ind in waters_generated: print ind
                                print "NOT_COPLANAR ELSE 1: waters_generated=", waters_generated
                    else:
                        #print "HYDROXYL/ETHER MODE"
                        lp_waters = self.Osp2(atom, master[0], master[1])
                        if self.verbose: print "lp_waters:"
                        for w in lp_waters:
                            if self.verbose:print w
                            waters_generated.append(w)
                        if self.verbose:
                            print "ELSE: waters_generated=", 
                            for ind in waters_generated: print ind
                            print "before adding else of not HYDROXYL, water_mates=",
                            for ind in water_mates: print ind
                    water_mates.append(waters_generated)
                    if self.verbose:
                        print "LAST 1: waters_generated=", 
                        for ind in waters_generated: print ind
                        print "LAST 1: water_mates =", 
                        for wm in water_mates : print wm
                        print "LAST 1: waters_generated=", waters_generated
                        print "LAST 1: water_mates =", water_mates 
                    numbering_stuff.append([int(position), len(waters_generated)-1])

            # NITROGEN #####################################################################################
            if atype == "NA":
                if len(master) == 1:
                    # calculate the Water vector
                    wet = self.hydro(master[0], atom) #3: @@REPAIR THIS APPARENT DUPLICATE??
                    waters_generated.append(wet)
                    water_mates.append(waters_generated)
                    if self.verbose:
                        print "NA 1: wet=", 
                        for ind in wet: print ind
                        print "NA 1: waters_generated=", 
                        for ind in waters_generated: print ind
                        print "NA 1: water_mates=", 
                        for ind in water_mates: print ind
                    numbering_stuff.append([int(position), len(waters_generated)-1])

                # nitrile mode
                if len(master) == 2:
                    wet = self.hydro(master[0], atom, master[1]) #4: @@REPAIR THIS APPARENT DUPLICATE4??
                    if self.verbose: print "nitrile mode: wet:", wet
                    #for ind in wet: print ind
                    waters_generated.append(wet)
                    if self.verbose: 
                        print "N: waters_generated:",
                        for ind in waters_generated: print ind
                    water_mates.append(waters_generated)
                    if self.verbose: 
                        print "N: water_mates:",
                        for ind in water_mates: print ind
                        print "NA 2: wet=", wet
                        print "NA 2: waters_generated=", waters_generated
                        print "NA 2: water_mates=", water_mates
                    numbering_stuff.append([int(position), len(waters_generated)-1])
                    if self.verbose:
                        print "nitrile numbering_stuff:",
                        for ind in numbering_stuff: print ind

                # tertiary amine HB receptor
                if len(master) == 3: 
                    #print "EXPERIMENTAL! potential tertiary amine"
                    master_center = self.mean3(master[0], master[1], master[2])
                    wet = self.hydro(master_center, atom) #5: @@REPAIR THIS APPARENT DUPLICATE5??
                    if self.verbose: 
                        print "tertiary amine: wet=",
                        for ind in wet: print ind
                    waters_generated.append(wet)
                    if self.verbose: 
                        print "tertiary amine: waters_generated=",
                        for ind in waters_generated: print ind
                    water_mates.append(waters_generated)
                    if self.verbose: 
                        print "tertiary amine: water_mates=",
                        for ind in water_mates: print ind
                        print "NA 3: wet=", wet
                        print "NA 3: waters_generated=", waters_generated
                        print "NA 3: water_mates=", water_mates
                    numbering_stuff.append([int(position), len(waters_generated)-1])
                    if self.verbose: print numbering_stuff[-1], ' added to numbering_stuff'

        ctr = 0
        for mates in water_mates:
            index = input.index(mates[0])
            if self.verbose: 
                print "mates=", mates
                print "mates[0]=", mates[0]
            line = ""
            for atom in mates:
                line += atom
            if self.verbose: print "@@ %d input line was = %s"%(index,  input[index])
            input[index] = line
            if self.verbose:
                print "now %d input line set to %s"%(index,  line)
                print "~~~numbering_stuff[%d]=%s"%(ctr, numbering_stuff[ctr])
            ctr+=1
            #if index==20: raise 'abc'
        count = 1
        # nicely split and format the lines that need to be splitted from
        # the previous addition of water molecules
        final = []
        if self.verbose: print "len(input)=", len(input), ' len(final)=', len(final)
        for line in input:
            line = line.split("\n")
            if self.verbose:
                print "ZZZZ: line=", 
                for  YYYY in line: print YYYY
            if len(line)==3 and line[-1]!='': 
                raise 'THREE!'
                line = line.remove(line[1])
            ct = 0
            for item in line:
                #@@@FIX THIS HACK!!!
                if not item == "" and item.find("WAT")<0 and item.find("W    99 1")<0:
                    final.append(item)
                    if len(item)>13 and item[13]=='W':
                        if self.verbose: print "after adding W ", item," len(final)=", len(final)
                    ct += 1
        if self.verbose: print "len(input)=", len(input), ' len(final)=', len(final)

        # renumbering the PDBQT atoms
        for line in final:
            if self.verbose: print "at top line=", line, " contains WAT ", line.find("WAT")>0
            if line.find("WAT")<0 and (line[0:4] == "ATOM" or line[0:6] == "HETATM"):
                value = "%4s" % count
                idx = final.index(line)
                final[idx] = line[0:7]+value+line[11:]
                count += 1
                if self.verbose: print "in loop line=", line

        # process the BRANCH numbers
        for line in final:
            idx = final.index(line)
            if "BRANCH" in line:
                #print "Processing line => ", line
                line = line.split()
                value1, value2 = int(line[1]), int(line[2])
                #print "BEFORE =>", value1, value2 
                addendum1 = 0
                addendum2 = 0
                for mark in numbering_stuff:
                    #print "Mark is :\t POSITION:", mark[0], " COUNT", mark[1]
                    if value1 > mark[0]:
                        addendum1 += mark[1]
                        #print "value 1 is bigger than mark[0]:", value1, " | ", mark[0]
                    if value2 > mark[0]:
                        addendum2 += mark[1]
                        #print "value 2 is bigger than mark[0]:", value2, " | ", mark[0]
                value1 += addendum1
                value2 += addendum2
                #print "AFTER  =>", value1, value2
                #print "- - - - - "
                final[idx] = line[0]+" "+str(value1)+" "+str(value2)

        # Writing the output
        if name is None:
            name = output+os.path.sep+name+"_HYDRO.pdbqt"

        try:
            hyd_ligand = open(name, 'w')
        except:
            print "# Error in saving the file %s #", name
            exit(1)

        count_waters = 0
        waters = {}
        for line in final:
        #        if not line == "":
            if line.split()[-1] == "W":
                count_waters += 1
                waters[line] = 1
            print >> hyd_ligand, line
        print " - %d waters added" % count_waters
        exit()    
                    
                

