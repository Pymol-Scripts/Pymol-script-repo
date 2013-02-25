'''
http://www.pymolwiki.org/index.php/inertia_tensor.py
 
(c) August 2010 by Mateusz Maciejewski
matt (at) mattmaciejewski . com

Licnese: MIT

'''

from pymol.cgo import *
from pymol import cmd


def tensor(selection, name="tensor", state=1, quiet=1):

    """
DESCRIPTION
    
    This script will draw the eigenvectors of the inertia tensor of the selection.

ARGUMENTS
    
    selection = string: selection for the atoms included in the tensor calculation

    name = string: name of the tensor object to be created {default: "tensor"}
    
    state = int: state/model in the molecule object used in the tensor calculation

EXAMPLE
    
    PyMOL> run inertia_tensor.py
    PyMOL> tensor molecule_object & i. 2-58+63-120 & n. n+ca+c, "tensor_model5_dom2", 5
    
NOTES
    
    Requires numpy.
    """


    import numpy


    def draw_axes(start,ends,radius=.2,name_obj="tensor"):
        radius = float(radius)
        size = radius*15.
        origin_offset = radius * -25.
        obj = [
            CYLINDER, start[0],start[1],start[2],ends[0][0]+start[0],ends[0][1]+start[1],ends[0][2]+start[2], radius, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,  # red Ixx
            CYLINDER, start[0],start[1],start[2],(-1)*ends[0][0]+start[0],(-1)*ends[0][1]+start[1],(-1)*ends[0][2]+start[2], radius, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
            CYLINDER, start[0],start[1],start[2],ends[1][0]+start[0],ends[1][1]+start[1],ends[1][2]+start[2], radius, 1.0, 1.0, 1.0, 0., 1.0, 0.,   # green Iyy
            CYLINDER, start[0],start[1],start[2],(-1)*ends[1][0]+start[0],(-1)*ends[1][1]+start[1],(-1)*ends[1][2]+start[2], radius, 1.0, 1.0, 1.0, 0., 1.0, 0.,
            CYLINDER, start[0],start[1],start[2],ends[2][0]+start[0],ends[2][1]+start[1],ends[2][2]+start[2], radius, 1.0, 1.0, 1.0, 0., 0.0, 1.0,  # blue Izz
            CYLINDER, start[0],start[1],start[2],(-1)*ends[2][0]+start[0],(-1)*ends[2][1]+start[1],(-1)*ends[2][2]+start[2], radius, 1.0, 1.0, 1.0, 0., 0.0, 1.0,

            ]

        cmd.load_cgo(obj,name_obj)



    totmass=0.0
    x_com,y_com,z_com=0,0,0
    I11,I12,I13,I21,I22,I23,I31,I32,I33=0,0,0,0,0,0,0,0,0

    model=cmd.get_model(selection, state)

    for a in model.atom:

            x_com+=a.coord[0]*a.get_mass()
            y_com+=a.coord[1]*a.get_mass()
            z_com+=a.coord[2]*a.get_mass()
            totmass+=a.get_mass()

    x_com /= totmass; y_com /= totmass; z_com /= totmass

    if not int(quiet):
        print
        print "Center of mass: "
        print
        print x_com, y_com, z_com

    I=[]

    for index in range(9):
        I.append(0)

    for a in model.atom:

            temp_x,temp_y,temp_z=a.coord[0],a.coord[1],a.coord[2]
            temp_x-=x_com; temp_y-=y_com; temp_z-=z_com

            I[0]+=a.get_mass()*(temp_y**2+temp_z**2)
            I[4]+=a.get_mass()*(temp_x**2+temp_z**2)
            I[8]+=a.get_mass()*(temp_x**2+temp_y**2)
            I[1]-=a.get_mass()*temp_x*temp_y
            I[3]-=a.get_mass()*temp_x*temp_y
            I[2]-=a.get_mass()*temp_x*temp_z
            I[6]-=a.get_mass()*temp_x*temp_z
            I[5]-=a.get_mass()*temp_y*temp_z
            I[7]-=a.get_mass()*temp_y*temp_z
   

    tensor = numpy.array([(I[0:3]),(I[3:6]),(I[6:9])])
    eigens = numpy.linalg.eig(tensor)
    vals,vects = numpy.linalg.eig(tensor) # they come out unsorted, so the command below is needed

    eig_ord = numpy.argsort(vals) # a thing to note is that here COLUMN i corrensponds to eigenvalue i.

    ord_vals = vals[eig_ord]
    ord_vects = vects[:,eig_ord].T

    if not int(quiet):
        print
        print "Inertia tensor x, y, z eigenvalues:"
        print
        print ord_vals
        print
        print "Inertia tensor x, y, z eigenvectors:"
        print
        print ord_vects


    start=[x_com,y_com,z_com]
    ends=[[10*ord_vects[0][0],10*ord_vects[0][1],10*ord_vects[0][2]],
          [10*ord_vects[1][0],10*ord_vects[1][1],10*ord_vects[1][2]],
          [10*ord_vects[2][0],10*ord_vects[2][1],10*ord_vects[2][2]]]

    draw_axes(start,ends,name_obj=name)


cmd.extend("tensor", tensor)
