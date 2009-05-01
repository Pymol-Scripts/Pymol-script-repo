from pymol import cmd
import pymol
import os
import re
import string
 
def __init__(self):
	cmd.extend('transform_odb', transform_odb)
 
# Creates a new object name from selection after transforming it with O-style matrix
# found in matrix_file
# Author: Mark Saper <saper@umich.edu>
 
def transform_odb( name, selection, matrix_file='matrix.odb',  transpose=0):
 
	# open the file for reading
	matrix_file = os.path.expanduser(matrix_file)
	matrix_file = os.path.expandvars(matrix_file)
	theInFile = open ( matrix_file,"r")
 
	# what we'll store the results in
	theMatrix = []
 
	# read every line in the file and ...
	for theCurrLine in theInFile.readlines():
	   if (theCurrLine) and (theCurrLine[0] != '!') and (theCurrLine[0] != '.'):
		  # if the line isn't blank, make a list of items seperated by tabs
		  theNewRow = string.split (theCurrLine)
		  # add it in the matrix
		  theMatrix.extend ( theNewRow )
 
	# change matrix to pymol unsupported format
 
	theMatrix = [ theMatrix[0], theMatrix[3], theMatrix[6], theMatrix[9],
					  theMatrix[1], theMatrix[4], theMatrix[7], theMatrix[10],
					  theMatrix [2], theMatrix [5], theMatrix[8], theMatrix[11], 
					  0.0, 0.0, 0.0, 0.0 ]
	theMatrix = [ float(x) for x in theMatrix]	
 
	# close the file
	theInFile.close ()
 
	r = cmd.create ( name, selection)
	r = cmd.transform_object( name, theMatrix, transpose=transpose)
 
	return r
 
cmd.extend('transform_odb', transform_odb)

