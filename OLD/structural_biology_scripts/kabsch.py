#!python
 
##############################################################################
#
# @SUMMARY: -- QKabsch.py.  A python implementation of the optimal superposition
#     of two sets of vectors as proposed by Kabsch 1976 & 1978.
#
# @AUTHOR: Jason Vertrees
# @COPYRIGHT: Jason Vertrees (C), 2005-2007
# @LICENSE: Released under GPL:
# This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA 
#
# DATE  : 2007-01-01
# REV   : 2
# REQUIREMENTS: numpy
#
#############################################################################
from array import *
 
# system stuff
import os
import copy
 
# pretty printing
import pprint
 
# for importing as a plugin into PyMol
from pymol import cmd
from pymol import stored
from pymol import selector
 
# using numpy for linear algebra
import numpy
 
def optAlign( sel1, sel2 ):
	"""
	optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
	Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
	Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA
 
	Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
	PyMol based upon your selections.
 
	By default, this program will optimally align the ALPHA CARBONS of the selections provided.
	To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
 
	@param sel1: First PyMol selection with N-atoms
	@param sel2: Second PyMol selection with N-atoms
	"""
	cmd.reset()
 
	# make the lists for holding coordinates
	# partial lists
	stored.sel1 = []
	stored.sel2 = []
	# full lists
	stored.mol1 = []
	stored.mol2 = []
 
	# -- CUT HERE
	sel1 = sel1 + " and N. CA"
	sel2 = sel2 + " and N. CA"
	# -- CUT HERE
 
	# Get the selected coordinates.  We
	# align these coords.
	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
 
	# get molecule name
	mol1 = cmd.identify(sel1,1)[0][0]
	mol2 = cmd.identify(sel2,1)[0][0]
 
	# Get all molecule coords.  We do this because
	# we have to rotate the whole molcule, not just
	# the aligned selection
	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
 
	# check for consistency
	assert( len(stored.sel1) == len(stored.sel2))
	L = len(stored.sel1)
	assert( L > 0 )
 
	# must alway center the two proteins to avoid
	# affine transformations.  Center the two proteins
	# to their selections.
	COM1 = numpy.sum(stored.sel1,axis=0) / float(L)
	COM2 = numpy.sum(stored.sel2,axis=0) / float(L)
	stored.sel1 = stored.sel1 - COM1
	stored.sel2 = stored.sel2 - COM2
 
	# Initial residual, see Kabsch.
	E0 = numpy.sum( numpy.sum(stored.sel1 * stored.sel1,axis=0),axis=0) + numpy.sum( numpy.sum(stored.sel2 * stored.sel2,axis=0),axis=0)
 
	#
	# This beautiful step provides the answer.  V and Wt are the orthonormal
	# bases that when multiplied by each other give us the rotation matrix, U.
	# S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
	V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(stored.sel2), stored.sel1))
 
	# we already have our solution, in the results from SVD.
	# we just need to check for reflections and then produce
	# the rotation.  V and Wt are orthonormal, so their det's
	# are +/-1.
	reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
 
	if reflect == -1.0:
		S[-1] = -S[-1]
		V[:,-1] = -V[:,-1]
 
	RMSD = E0 - (2.0 * sum(S))
	RMSD = numpy.sqrt(abs(RMSD / L))
 
	#U is simply V*Wt
	U = numpy.dot(V, Wt)
 
	# rotate and translate the molecule
	stored.sel2 = numpy.dot((stored.mol2 - COM2), U)
	stored.sel2 = stored.sel2.tolist()
	# center the molecule
	stored.sel1 = stored.mol1 - COM1
	stored.sel1 = stored.sel1.tolist()
 
	# let PyMol know about the changes to the coordinates
	cmd.alter_state(1,mol1,"(x,y,z)=stored.sel1.pop(0)")
	cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
 
	print "RMSD=%f" % RMSD
 
	# make the alignment OBVIOUS
	cmd.hide('everything')
	cmd.show('ribbon', sel1 + ' or ' + sel2)
	cmd.color('gray70', mol1 )
	cmd.color('paleyellow', mol2 )
	cmd.color('red', 'visible')
	cmd.show('ribbon', 'not visible')
	cmd.center('visible')
	cmd.orient()
	cmd.zoom('visible')
 
cmd.extend("optAlign", optAlign)
