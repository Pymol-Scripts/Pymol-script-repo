
'''
Calculate and display the dimensions of a protein.


This is a first version, please use at your own risk!

REQUIREMENTS

numpy (http://numpy.scipy.org) that should be built into the newers versions of Pymol

(c) Pablo Guardado Calvo

Based on "inertia_tensor.py" (c) 2010 by Mateusz Maciejewski

License: MIT

'''

from __future__ import print_function

__author__  = 'Pablo Guardado Calvo'
__version__ = '0.2'
__email__   = 'pablo.guardado (at) gmail.com'
__date__    = '13/08/2015'
__modification_date__ = '05/02/2018'
__modification_reason__ = 'Error in the code produced sometimes inverted structures'

###########################################################################################################################################################
# USAGE
#
# The idea behing this script is to calculate an aproximate minimal bounding box to extract the cell dimensions of a protein. To calculate the minimal bounding
# is not trivial and usually the Axis Aligned Bounding Box (AABB) does not show up the real dimensions of the protein. This script calculates the inertia tensor
# of the object, extract the eigenvalues and use them to rotate the molecule (using as rotation matrix the transpose of the eigenvalues matrix). The result is that
# the molecule is oriented with the inertia axis aligned with the cartesian axis. A new Bounding Box is calculated that is called Inertia Axis Aligned Bounding Box
#(IABB), whose volume is always lower than AABB volume, and in many cases will correspond with the lowest volume. Of course, maybe it exists another Bounding Box
# with a lower volume (the minimal Bounding Box).
#
# As always with these type of things, you have to use at your own risk. I did not try all the possible combinations, but if you find a bug, do
# not hesitate to contact me (pablo.guardado (at) gmail.com) or try to modify the code for yourself to correct it.
#
# To load the script just type:
#
# run path-to-the-script/Draw_Protein_Dimensions.py
#
# or if you want something more permanent add the previous line to your .pymolrc file
#
# The script works just typing:
#
# draw_Protein_Dimensions selection
#
# This will draw the cell dimensions of your selection based on a IABB. It also generates the IABB box and the inertia axis, you just need to do "show cgo" to display them.
#
# You could also try:
#
# draw_BB selection
#
# This will draw the AABB and IABB boxes with their cell dimensions and show in the command line their volumes, you can compare both of them.
############################################################################################################################################################

from pymol import cmd, cgo
from pymol.cgo import *
import numpy
from random import randint

def matriz_inercia(selection):
	'''
	DESCRIPTION

	The method calculates the mass center, the inertia tensor and the eigenvalues and eigenvectors
	for a given selection. Mostly taken from inertia_tensor.py
	'''

	model = cmd.get_model(selection)
	totmass = 0.0
	x,y,z = 0,0,0
	for a in model.atom:
		m = a.get_mass()
		x += a.coord[0]*m
		y += a.coord[1]*m
		z += a.coord[2]*m
		totmass += m
	global cM
	cM = numpy.array([x/totmass, y/totmass, z/totmass])


	I = []
	for index in range(9):
		I.append(0)

	for a in model.atom:
		temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
		temp_x -= x
		temp_y -= y
		temp_z -= z

		I[0] += a.get_mass() * (temp_y**2 + temp_z**2)
		I[1] -= a.get_mass() * temp_x * temp_y
		I[2] -= a.get_mass() * temp_x * temp_z
		I[3] -= a.get_mass() * temp_x * temp_y
		I[4] += a.get_mass() * (temp_x**2 + temp_z**2)
		I[5] -= a.get_mass() * temp_y * temp_z
		I[6] -= a.get_mass() * temp_x * temp_z
		I[7] -= a.get_mass() * temp_y * temp_z
		I[8] += a.get_mass() * (temp_x**2 + temp_y**2)

	global tensor
	tensor = numpy.array([(I[0:3]), (I[3:6]), (I[6:9])])

	global autoval, autovect, ord_autoval, ord_autovect
	autoval, autovect = numpy.linalg.eig(tensor)
	auto_ord = numpy.argsort(autoval)
	ord_autoval = autoval[auto_ord]
	ord_autovect_complete = autovect[:, auto_ord].T
	ord_autovect = numpy.around(ord_autovect_complete, 3)

	return ord_autoval


def draw_inertia_axis(selection):
	'''
	DESCRIPTION

	This method draw the inertia axis calculated with the method matriz_inercia.
	'''

	matriz_inercia(selection)
	axis1 = ord_autovect[0]
	x1, y1, z1 = cM[0], cM[1], cM[2]
	x2, y2, z2 = cM[0]+50*axis1[0], cM[1]+50*axis1[1], cM[2]+50*axis1[2]
	eje1 = [cgo.CYLINDER, x1, y1, z1, x2, y2, z2, 0.6, 1, 0, 0, 1, 0, 0, 0.0]
	cmd.load_cgo(eje1, 'Inertia_Axis1')
	axis2 = ord_autovect[1]
	x3, y3, z3 = cM[0]+40*axis2[0], cM[1]+40*axis2[1], cM[2]+40*axis2[2]
	eje1 = [cgo.CYLINDER, x1, y1, z1, x3, y3, z3, 0.6, 1, 0.5, 0, 1, 0.5, 0, 0.0]
	cmd.load_cgo(eje1, 'Inertia_Axis2')
	axis4 = ord_autovect[2]
	x4, y4, z4 = cM[0]+30*axis4[0], cM[1]+30*axis4[1], cM[2]+30*axis4[2]
	eje1 = [cgo.CYLINDER, x1, y1, z1, x4, y4, z4, 0.6, 1, 1, 0, 1, 1, 0, 0.0]
	cmd.load_cgo(eje1, 'Inertia_Axis3')

def translacion_cM(selection):
	'''
	DESCRIPTION

	Translate the center of mass of the molecule to the origin.
	'''
	model = cmd.get_model(selection)
	totmass = 0.0
	x,y,z = 0,0,0
	for a in model.atom:
		m = a.get_mass()
		x += a.coord[0]*m
		y += a.coord[1]*m
		z += a.coord[2]*m
		totmass += m
	cM = numpy.array([x/totmass, y/totmass, z/totmass])
	trans_array = ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -cM[0], -cM[1], -cM[2], 1])
	model_trans = cmd.transform_selection(selection, trans_array)


def rotacion_orig(selection):
	'''
	DESCRIPTION

	Find the proper rotation matrix, i.e. the transpose of the matrix formed by the eigenvectors of the inertia tensor
	'''

	translacion_cM(selection)
	matriz_inercia(selection)
	global transf, transf_array, ord_autovect_array, transf_array_print
	ord_autovect_array = numpy.array([[ord_autovect[0][0], ord_autovect[0][1], ord_autovect[0][2]],
	                                  [ord_autovect[1][0], ord_autovect[1][1], ord_autovect[1][2]],
					  [ord_autovect[2][0], ord_autovect[2][1], ord_autovect[2][2]]])
	if numpy.linalg.det(ord_autovect_array) < 0.:
		ord_autovect_array = numpy.array([[ord_autovect[2][0], ord_autovect[2][1], ord_autovect[2][2]],
	                                  	  [ord_autovect[1][0], ord_autovect[1][1], ord_autovect[1][2]],
	 				  	  [ord_autovect[0][0], ord_autovect[0][1], ord_autovect[0][2]]])
	transf = numpy.transpose(ord_autovect_array)
	transf_array = numpy.array([transf[0][0], transf[0][1], transf[0][2], 0,
	                            transf[1][0], transf[1][1], transf[1][2], 0,
				    transf[2][0], transf[2][1], transf[2][2], 0,
				    0, 0, 0, 1])


def transformar(selection):
	'''
	DESCRIPTION

	Rotate the molecule and draw the inertia axis.
	'''

	rotacion_orig(selection)
	model_rot = cmd.transform_selection(selection, transf_array, homogenous=0, transpose=1);
	draw_inertia_axis(selection)



def draw_AABB(selection):
        """
        DESCRIPTION
        For a given selection, draw the Axes Aligned bounding box around it without padding. Code taken and modified from DrawBoundingBox.py.

        """

 	AA_original = selection + "_original"
	model_orig = cmd.create(AA_original, selection)

        ([min_X, min_Y, min_Z],[max_X, max_Y, max_Z]) = cmd.get_extent(AA_original)

	print("The Axis Aligned Bounding Box (AABB) dimensions are (%.2f, %.2f, %.2f)" % (max_X-min_X, max_Y-min_Y, max_Z-min_Z))
	print("The Axis Aligned Bounding Box (AABB) volume is %.2f A3" % ((max_X-min_X)*(max_Y-min_Y)*(max_Z-min_Z)))


        min_X = min_X
        min_Y = min_Y
        min_Z = min_Z
        max_X = max_X
        max_Y = max_Y
        max_Z = max_Z

        boundingBox = [
                LINEWIDTH, float(2),

                BEGIN, LINES,
                COLOR, float(1), float(1), float(0),

                VERTEX, min_X, min_Y, min_Z,
                VERTEX, min_X, min_Y, max_Z,
                VERTEX, min_X, max_Y, min_Z,
                VERTEX, min_X, max_Y, max_Z,
                VERTEX, max_X, min_Y, min_Z,
                VERTEX, max_X, min_Y, max_Z,
                VERTEX, max_X, max_Y, min_Z,
                VERTEX, max_X, max_Y, max_Z,
                VERTEX, min_X, min_Y, min_Z,
                VERTEX, max_X, min_Y, min_Z,
                VERTEX, min_X, max_Y, min_Z,
                VERTEX, max_X, max_Y, min_Z,
                VERTEX, min_X, max_Y, max_Z,
                VERTEX, max_X, max_Y, max_Z,
                VERTEX, min_X, min_Y, max_Z,
                VERTEX, max_X, min_Y, max_Z,
                VERTEX, min_X, min_Y, min_Z,
                VERTEX, min_X, max_Y, min_Z,
                VERTEX, max_X, min_Y, min_Z,
                VERTEX, max_X, max_Y, min_Z,
                VERTEX, min_X, min_Y, max_Z,
                VERTEX, min_X, max_Y, max_Z,
                VERTEX, max_X, min_Y, max_Z,
                VERTEX, max_X, max_Y, max_Z,

		END
	]

	p0 = '_0'  + str(randint(0, 100))
	p1 = '_1' + str(randint(0, 100))
	p2 = '_2' + str(randint(0, 100))
	p3 = '_3' + str(randint(0, 100))
	cmd.pseudoatom (pos=[min_X, min_Y, min_Z], object=p0)
	cmd.pseudoatom (pos=[min_X, min_Y, max_Z], object=p1)
	cmd.pseudoatom (pos=[min_X, max_Y, min_Z], object=p2)
	cmd.pseudoatom (pos=[max_X, min_Y, min_Z], object=p3)
	cmd.distance(None, p0, p3)
	cmd.distance(None, p0, p2)
	cmd.distance(None, p0, p1)
	cmd.hide("nonbonded")

	boxName = "box_AABB_"  + str(randint(0, 100))
        cmd.load_cgo(boundingBox,boxName)
        return boxName


def draw_IABB(selection):
        """
        DESCRIPTION
        For a given selection, draw the Inertia Axes Aligned bounding box around it without padding. Code taken and modified from DrawBoundingBox.py.

        """

	transformar(selection)

        ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

	print("The Inertia Axis Aligned Bounding Box (IABB) dimensions are (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))
	print("The Inertia Axis Aligned Bounding Box (IABB) volume is %.2f A3" % ((maxX-minX)*(maxY-minY)*(maxZ-minZ)))


        minX = minX
        minY = minY
        minZ = minZ
        maxX = maxX
        maxY = maxY
        maxZ = maxZ

        boundingBox = [
                LINEWIDTH, float(2),

                BEGIN, LINES,
                COLOR, float(1), float(0), float(0),

                VERTEX, minX, minY, minZ,
                VERTEX, minX, minY, maxZ,
                VERTEX, minX, maxY, minZ,
                VERTEX, minX, maxY, maxZ,
                VERTEX, maxX, minY, minZ,
                VERTEX, maxX, minY, maxZ,
                VERTEX, maxX, maxY, minZ,
                VERTEX, maxX, maxY, maxZ,
                VERTEX, minX, minY, minZ,
                VERTEX, maxX, minY, minZ,
                VERTEX, minX, maxY, minZ,
                VERTEX, maxX, maxY, minZ,
                VERTEX, minX, maxY, maxZ,
                VERTEX, maxX, maxY, maxZ,
                VERTEX, minX, minY, maxZ,
                VERTEX, maxX, minY, maxZ,
                VERTEX, minX, minY, minZ,
                VERTEX, minX, maxY, minZ,
                VERTEX, maxX, minY, minZ,
                VERTEX, maxX, maxY, minZ,
                VERTEX, minX, minY, maxZ,
                VERTEX, minX, maxY, maxZ,
                VERTEX, maxX, minY, maxZ,
                VERTEX, maxX, maxY, maxZ,

                END
        ]

 	p4 = '_4' + str(randint(0, 100))
	p5 = '_5' + str(randint(0, 100))
	p6 = '_6' + str(randint(0, 100))
	p7 = '_7' + str(randint(0, 100))
	cmd.pseudoatom (pos=[minX, minY, minZ], object=p4)
	cmd.pseudoatom (pos=[minX, minY, maxZ], object=p5)
	cmd.pseudoatom (pos=[minX, maxY, minZ], object=p6)
	cmd.pseudoatom (pos=[maxX, minY, minZ], object=p7)
	cmd.distance(None, p4, p7)
	cmd.distance(None, p4, p6)
	cmd.distance(None, p4, p5)
	cmd.hide("nonbonded")

	boxName = "box_IABB_" + str(randint(0, 100))
        cmd.load_cgo(boundingBox,boxName)
        return boxName


def draw_BB(selection):
	draw_AABB(selection)
	draw_IABB(selection)

def draw_Protein_Dimensions(selection):
	draw_IABB(selection)
	cmd.hide("cgo")

cmd.extend ("draw_Protein_Dimensions", draw_Protein_Dimensions)
cmd.extend ("draw_BB", draw_BB)







