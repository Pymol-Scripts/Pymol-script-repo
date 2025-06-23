'''
Calculate and display the rotation axis and the translation vector for a given transformation.

http://pymolwiki.org/index.php/RotationAxis

This is a first version, please use at your own risk!

REQUIREMENTS

numpy (http://numpy.scipy.org) that should be built into the newers versions of Pymol

'''

from __future__ import print_function

__author__  = 'Pablo Guardado Calvo'
__version__ = '0.1'
__email__   = 'pablo.guardado (at) gmail.com'

###########################################################################################################################################################
# USAGE
#
# The idea behind this script is to align two molecules/domains/chains/selections (using cmd.super) and extract the trasformation (TTT) matrix (T).
# From this matrix we obtain the direction of the rotation axis, as well as a point and we create a cgo object representing the axis. The script
# generates two measures: one in the graphical screen (the rotation axis value, and the norm of the tranlation vector along the rotation axis) and
# some basic information in the command-line (the transformation matrix, the rotation angle, distance between centroids, some pml lines you can use
# to reproduce the axis...)
#
# As always with these type of things, you have to use at your own risk. I did not try all the possible combination, but if you find a bug, do
# not hesitate to contact me (pablo.guardado (at) gmail.com) or try to modify the code for yourself to correct it.
#
# To load the script just type:
#
# run path-to-the-script/draw_rotation_axis.py
#
# or if you want something more permanent add the previous line to your .pymolrc file
#
# The script works just typing:
#
# draw_axis('selection1', 'selection2')
#
# Please, pay attention to the apostrophes around the selections, you MUST use them. Also works with chains:
#
# draw_axis('chain A', 'chain B')
#
# Also, you can play a bit with the lenght, width and colour of the axis you are going to generate.
#
# draw_axis('selection1', 'selection2', scale_factor, width, r1, g1, b1, r2, g2, b2)
#
# scale_factor = to control the lenght of the axis, the default is 20
# width = to control the width of the axis. Default is 0.6
# r1, g1, b1 = first RGB colour code. Default is 1, 1, 1
# r2, g2, b2 = second RGB colour code to create the gradient. Default is 1, 0, 0.
# To create a single colour axis, just made r1,g1,b1=r2,g2,b2
#
############################################################################################################################################################


from pymol import cmd, cgo
import math
import numpy


def transf_matrix(chA, chB):
    '''
    DESCRIPTION

    Align two selections/chains, and returns the transformation matrix. I used super to carry out the alignment, likely is possible to use cmd.align and
    is going to be a bit faster, but I think is not going to work well with low-sequence-identity alignments.

    '''
    cmd.create('working', chA)
    cmd.super('working', chB)
    T = cmd.get_object_matrix('working')
    global cmW
    cmW = center_of_Mass('working')
    cmd.delete('working')
    return T


def center_of_Mass(selection):
    '''
    DESCRIPTION

    Calculates the center of mass of a given selection

    '''
    model= cmd.get_model(selection)
    x,y,z=0,0,0
    totmass = 0
    for a in model.atom:
        m = a.get_mass()
        x+= a.coord[0]*m
        y+= a.coord[1]*m
        z+= a.coord[2]*m
        totmass += m
    cM = numpy.array([x/totmass, y/totmass, z/totmass])
    return cM

def direction_cosines(chA, chB):
    '''
    DESCRIPTION

    Calculates the direction cosines of the rotation axis from the transformation matrix.

    '''
    t=transf_matrix(chA, chB)
    a1= (t[6]-t[9])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    b1= (t[8]-t[2])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    c1= (t[1]-t[4])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    axis = numpy.array([a1, b1, c1])
    return axis

def angle_axis(chA, chB):
    '''
    DESCRIPTION

    Calculates the rotation angle from the transformation matrix

    '''
    t=transf_matrix(chA, chB)
    angle_rad = math.acos((t[0]+t[5]+t[10]-1)/2)
    return angle_rad

def center_of_rot(chA, chB):
    '''
    DESCRIPTION

    Calculates the center of rotation of the axis

    '''
    cm_Working=center_of_Mass(chA)
    cm_Reference=cmW
    u=direction_cosines(chA, chB)[0]
    u2=u**2
    v=direction_cosines(chA, chB)[1]
    v2=v**2
    w=direction_cosines(chA, chB)[2]
    w2=w**2
    cos_theta=numpy.cos(angle_axis(chA, chB))
    sin_theta=numpy.sin(angle_axis(chA, chB))
    fx=cm_Working[0]
    fy=cm_Working[1]
    fz=cm_Working[2]
    x=cm_Reference[0]
    y=cm_Reference[1]
    z=cm_Reference[2]

    T11 = (v2 + w2)*(1-cos_theta)
    T12 = (w*sin_theta)-((u*v)*(1-cos_theta))
    T13 = -(v*sin_theta)-(u*w*(1-cos_theta))
    T14 =fx-((((u2*x)+((u*v)*y)+((u*w)*z))*(1-cos_theta))+(x*cos_theta)+((-(w*y)+(v*z))*sin_theta))
    T21 = -(w*sin_theta)-((u*v)*(1-cos_theta))
    T22 = (u2 + w2)*(1-cos_theta)
    T23 =  (u*sin_theta)-(w*v*(1-cos_theta))
    T24 =fy-((((v*u*x)+(v2*y)+(v*w*z))*(1-cos_theta))+(y*cos_theta)+(((w*x)-(u*z))*sin_theta))
    T31 = (v*sin_theta)-(w*u*(1-cos_theta))
    T32 = -(u*sin_theta)-(w*v*(1-cos_theta))
    T33 = (u2 + v2)*(1-cos_theta)
    T34 =fz-(((((u*x)*w)+((v*y)*w)+(w2*z))*(1-cos_theta))+(z*cos_theta)+((-(v*x)+(u*y))*sin_theta))

    term_lig = numpy.array([[T11, T12, T13], [T21, T22, T23], [T31, T32, T33]])
    term_ind = numpy.array([T14, T24, T34])

    sol_lstsq = numpy.linalg.lstsq(term_lig, term_ind)
    sol = sol_lstsq[0]

    return sol

def nearest_point_to_axis(chA, chB):
    '''
    DESCRIPTION

    Calculates the nearest point of the axis, I use it to create the cgo object.

    '''
    cmA=center_of_Mass(chA)
    cmB=cmW
    cmAver=(cmB+cmA)/2
    vector=numpy.array([(cmB[0]-cmA[0]), (cmB[1]-cmA[1]), (cmB[2]-cmA[2])])
    moduli_vector=numpy.linalg.norm(vector)
    vector_director=numpy.array([(cmB[0]-cmA[0])/moduli_vector, (cmB[1]-cmA[1])/moduli_vector, (cmB[2]-cmA[2])/moduli_vector])
    axis1= direction_cosines(chA, chB)
    sol=center_of_rot(chA, chB)
    term_lig2=numpy.array([[vector_director[0], vector_director[1], vector_director[2], 0], [1, 0, 0, -axis1[0]], [0, 1, 0, -axis1[1]], [0, 0, 1, -axis1[2]]])
    term_ind2=numpy.array([(cmAver[0]*(vector_director[0]))+(cmAver[1]*(vector_director[1]))+(cmAver[2]*(vector_director[2])), sol[0], sol[1], sol[2]])
    term_j=(cmAver[0]*vector_director[0])+(cmAver[1]*vector_director[1])+(cmAver[2]*vector_director[2])
    suma_vect_director=vector_director+axis1
    term_ji=(cmAver[0]*suma_vect_director[0])+(cmAver[1]*suma_vect_director[1])+(cmAver[2]*suma_vect_director[2])
    if numpy.dot(vector_director, axis1) != 0:
        t = ((-numpy.dot(vector_director, sol))+term_j)/numpy.dot(vector_director, axis1)
    else:
        t = ((-numpy.dot(suma_vect_director, sol))+term_ji)/numpy.dot(suma_vect_director, axis1)
    p =  [sol[0]+axis1[0]*t, sol[1]+axis1[1]*t, sol[2]+axis1[2]*t]

    return p

def proyeccion_centroide(selection, chA, chB):
    '''
    DESCRIPTION

    Calculates the proyection of the mass center for the working molecule before being aligned with the reference molecule. For representation purpuses.

    '''
    axis1=numpy.array([direction_cosines(chA, chB)[0], direction_cosines(chA, chB)[1], direction_cosines(chA, chB)[2]])
    sol=center_of_rot(chA, chB)
    cmSel=center_of_Mass(selection)
    t_cen=numpy.dot(cmSel, axis1)-numpy.dot(sol, axis1)
    proy_cen= [sol[0]+(t_cen*axis1[0]), sol[1]+(t_cen*axis1[1]), sol[2]+(t_cen*axis1[2])]
    return proy_cen

def proyeccion_centroide_working(chA, chB):
    '''
    DESCRIPTION

    Calculates the proyection of the mass center for working molecule after being aligned with the reference molecule. For representation purpuses.

    '''
    axis1=numpy.array([direction_cosines(chA, chB)[0], direction_cosines(chA, chB)[1], direction_cosines(chA, chB)[2]])
    sol=center_of_rot(chA, chB)
    cmSel=cmW
    t_cen=numpy.dot(cmSel, axis1)-numpy.dot(sol, axis1)
    proy_cen= [sol[0]+(t_cen*axis1[0]), sol[1]+(t_cen*axis1[1]), sol[2]+(t_cen*axis1[2])]
    return proy_cen

def print_information(T, axis1, angle_degrees,  moduli_vector, obj, x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2, modu_tr=0):
    '''
    DESCRIPTION

    Print to basic information to the screen.
    '''
    print("#################################################################################################")
    print("Transformation (TTT) matrix")
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[0], T[1], T[2], T[3]))
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[4], T[5], T[6], T[7]))
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[8], T[9], T[10], T[11]))
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[12], T[13], T[14], T[15]))
    print(".................................................................................................")
    print("")
    print("The direction cosines of the rotation axis is: %3.2f, %3.2f, %3.2f" %(axis1[0], axis1[1], axis1[2]))
    print("The angle of rotation is %3.2f degrees" %(angle_degrees))
    print("The lenght of the translation vector along the rotation axis is %3.2f Angstroms" %(modu_tr))
    print("The distance between mass centers is %3.2f Angstroms" %(moduli_vector))
    print(".................................................................................................")
    print("")
    print("Lines to be used in a pml script to generate the axis")
    print("")
    print("CYLINDER, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, 0.0" %(x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2))
    print("cmd.load_cgo(obj, %3.2f)" %(angle_degrees))
    print("")
    print("#################################################################################################")

def draw_axis(chA, chB, scale_factor=20, w=0.6, r1=1, g1=1, b1=1, r2=1, g2=0, b2=0):
    T = transf_matrix(chA, chB)
    angle=angle_axis(chA, chB)
    angle_degrees=(angle*180)/math.pi
    axis1=[direction_cosines(chA, chB)[0], direction_cosines(chA, chB)[1], direction_cosines(chA, chB)[2]]
    p = nearest_point_to_axis(chA, chB)
    x1, y1, z1 = p[0] + (3*scale_factor*axis1[0]), p[1] + (3*scale_factor*axis1[1]), p[2] + (3*scale_factor*axis1[2])
    x2, y2, z2 = p[0] - (3*scale_factor*axis1[0]), p[1] - (3*scale_factor*axis1[1]), p[2] - (3*scale_factor*axis1[2])
    obj = [cgo.CYLINDER, x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2, 0.0]
    cmd.load_cgo(obj, angle_degrees)

    cmA=center_of_Mass(chA)
    cmB=cmW
    cmAver=(cmB+cmA)/2
    vector=numpy.array([(cmB[0]-cmA[0]), (cmB[1]-cmA[1]), (cmB[2]-cmA[2])])
    moduli_vector=numpy.linalg.norm(vector)
    vector_director=numpy.array([(cmB[0]-cmA[0])/moduli_vector, (cmB[1]-cmA[1])/moduli_vector, (cmB[2]-cmA[2])/moduli_vector])
    pC_A = proyeccion_centroide(chA, chA, chB)
    pC_B = proyeccion_centroide_working(chA, chB)


    trans_vector = numpy.array([(pC_B[0]-pC_A[0]), (pC_B[1]-pC_A[1]), (pC_B[2]-pC_A[2])])
    modu_tr = numpy.linalg.norm(trans_vector)
    rota_centroid_rad=numpy.dot(vector_director, axis1)
    rota_centroid = (rota_centroid_rad*180)/math.pi
    rota_centroid_absol_0= numpy.absolute(rota_centroid)
    rota_centroid_absol=round(rota_centroid_absol_0,2)


    if rota_centroid_absol == 0.00:
        p1 = '_1'
        p2 = '_2'
        p3 = '_3'
        cmd.pseudoatom (pos=[cmA[0], cmA[1], cmA[2]], object=p1)
        cmd.pseudoatom (pos=[pC_A[0], pC_A[1], pC_A[2]], object=p2)
        cmd.pseudoatom (pos=[cmB[0], cmB[1], cmB[2]], object=p3)
        cmd.angle(None, p1, p2, p3)
        print_information(T, axis1, angle_degrees,  moduli_vector, obj, x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2)

    if rota_centroid_absol != 0:
        p1 = '_1'
        p2 = '_2'
        p3 = '_3'
        p4 = '_4'
        cmd.pseudoatom (pos=[cmA[0], cmA[1], cmA[2]], object=p1)
        cmd.pseudoatom (pos=[pC_A[0], pC_A[1], pC_A[2]], object=p2)
        cmd.pseudoatom (pos=[pC_B[0], pC_B[1], pC_B[2]], object=p3)
        cmd.pseudoatom (pos=[cmB[0], cmB[1], cmB[2]], object=p4)
        cmd.dihedral(None, p1, p2, p3, p4)
        cmd.distance(None, p2, p3)
        print_information(T, axis1, angle_degrees,  moduli_vector, obj, x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2, modu_tr)

    cmd.create('working', chA)
    cmd.super('working', chB)
