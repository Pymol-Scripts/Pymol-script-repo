"""
    Parser for TOPOLOGY.xml

    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Copyright (c) 2002-2011, Jens Erik Nielsen, University College Dublin; 
    Nathan A. Baker, Battelle Memorial Institute, Developed at the Pacific 
    Northwest National Laboratory, operated by Battelle Memorial Institute, 
    Pacific Northwest Division for the U.S. Department Energy.; 
    Paul Czodrowski & Gerhard Klebe, University of Marburg.

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
        * Neither the names of University College Dublin, Battelle Memorial Institute,
          Pacific Northwest National Laboratory, US Department of Energy, or University
          of Marburg nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

    ----------------------------

"""

__date__ = "12 November 2008"
__author__ = "Nathan Baker, Yong Huang"


TOPOLOGYPATH = "TOPOLOGY.xml"

from sys import stderr
from xml import sax

class TopologyHandler(sax.ContentHandler):
	""" Handler for XML-based topology files.  Assumes the following hierarchy of tags:
	topology
	-->residue
	   |-->reference
	   |-->titrationstate
	       |-->tautomer
	           |-->conformer
	"""
	def __init__(self):
		self.currentElement = None
		self.currentAtom = None
		self.currentDihedral = None
		self.currentReference = None
		self.currentResidue = None
		self.currentTitrationState = None
		self.currentTautomer = None
		self.currentConformer = None
		self.currentConformerAdd = None
		self.currentConformerRemove = None
		self.residues = []
		self.incomplete = 0
		
	def startElement(self, tagName, attributes):
		if not self.incomplete: 
			#print "Processing %s start tag" % tagName
			if tagName == "topology":
				pass
			elif tagName == "residue":
				if self.currentResidue != None:
					print "** Overwriting current TopologyResidue object!"
				self.currentResidue = TopologyResidue(self)
			elif tagName == "reference":
				if self.currentReference != None:
					print "** Overwriting current TopologyReference object!"
				self.currentReference = TopologyReference(self.currentResidue)
			elif tagName == "titrationstate":
				if self.currentTitrationState != None:
					print "** Overwriting current TopologyTitrationState object!"
				self.currentTitrationState = TopologyTitrationState(self.currentResidue)
			elif tagName == "tautomer":
				if self.currentTautomer != None:
					print "** Overwriting current Tautomer object!"
				self.currentTautomer = TopologyTautomer(self.currentTitrationState)
			elif tagName == "conformer":
				if self.currentConformer != None:
					print "** Overwriting current Conformer object!"
				self.currentConformer = TopologyConformer(self.currentTautomer)
			elif tagName == "name":
				self.currentElement = tagName
			elif tagName == "atom":
				if self.currentConformerAdd != None:
					#print "    Adding atom to conformerAdd..."
					self.currentAtom = TopologyAtom(self.currentConformerAdd)
				elif self.currentConformerRemove != None:
					#print "    Adding atom to conformerRemove..."
					self.currentAtom = TopologyAtom(self.currentConformerRemove)
				elif self.currentReference != None:
					#print "    Adding atom to reference..."
					self.currentAtom = TopologyAtom(self.currentReference)
				else:
					print "** Don't know what to do with this atom!"
			elif tagName == "x":
				self.currentElement = tagName
			elif tagName == "y":
				self.currentElement = tagName
			elif tagName == "z":
				self.currentElement = tagName
			elif tagName == "bond":
				self.currentElement = tagName
			elif tagName == "altname":
				self.currentElement = tagName
			elif tagName == "dihedral":
				self.currentElement = tagName
				if self.currentConformerAdd != None:
					#print "    Adding dihedral to conformerAdd..."
					self.currentDihedral = TopologyDihedral(self.currentConformerAdd)
				elif self.currentConformerRemove != None:
					#print "    Adding dihedral to conformerRemove..."
					self.currentDihedral = TopologyDihedral(self.currentConformerRemove)
				elif self.currentReference != None:
					#print "    Adding dihedral to reference..."
					self.currentDihedral = TopologyDihedral(self.currentReference)
				else:
					print "** Don't know what to do with this dihedral!"
			elif tagName == "add":
				self.currentConformerAdd = TopologyConformerAdd(self.currentConformer)
			elif tagName == "remove":
				#print "currentConformer: %s" % (self.currentConformer)
				self.currentConformerRemove = TopologyConformerRemove(self.currentConformer)
			elif tagName == "incomplete":
				#print "incomplete state encounted, skipping!"
				self.incomplete = 1
			else:
				print "** NOT handling %s start tag" % tagName
			
	def endElement(self, tagName):
		if not self.incomplete:
			#print "Processing %s end tag" % tagName
			self.currentElement == None
			if tagName == "x":
				pass
			elif tagName == "y":
				pass
			elif tagName == "z":
				pass
			elif tagName == "name":
				pass
			elif tagName == "bond":
				pass
			elif tagName == "altname":
				pass
			elif tagName == "atom":
				self.currentAtom = None
			elif tagName == "dihedral":
				self.currentDihedral = None
			elif tagName == "reference":
				self.currentReference = None
			elif tagName == "add":
				self.currentConformerAdd = None
			elif tagName == "remove":
				self.currentConformerRemove = None
			elif tagName == "titrationstate":
				residue = self.currentTitrationState.topologyResidue
				#print "Residue %s has titration states:  %s" % (residue.name, residue.titrationStates)
				self.currentTitrationState = None
			elif tagName == "conformer":
				tautomer = self.currentConformer.topologyTautomer
				#print "Tautomer %s has conformers:  %s" % (tautomer.name, tautomer.conformers)
				self.currentConformer = None
			elif tagName == "tautomer":
				titrationState = self.currentTautomer.topologyTitrationState
				#print "Titration state %s has tautomers:  %s" % (titrationState.name, titrationState.tautomers)
				self.currentTautomer = None
			elif tagName == "residue":
				self.currentResidue = None
			elif tagName == "topology":
				pass
			else:
				print "** NOT handling %s end tag" % tagName
		else:
			if tagName == "incomplete":
				self.incomplete = 0

			
	def characters(self, text):
		if text.isspace():  return
		
		if not self.incomplete:
			if self.currentElement == "name":
				# Processing a name based on current context
				if self.currentAtom != None:
					#print "    Setting atom name to %s" % text
					self.currentAtom.name = text
				elif self.currentConformer != None:
					#print "    Setting conformer name to %s" % text
					self.currentConformer.name = text
				elif self.currentTautomer != None:
					#print "    Setting tautomer name to %s" % text
					self.currentTautomer.name = text
				elif self.currentTitrationState != None:
					#print "    Setting titration state name to %s" % text
					self.currentTitrationState.name = text
				elif self.currentResidue != None:
					#print "    Setting residue name to %s" % text
					self.currentResidue.name = text
				else:
					print "    *** Don't know what to do with name %s!" % text
			elif self.currentElement == "x":
				#print "    Setting atom x coordinate to %s" % text
				self.currentAtom.x = float(text)
			elif self.currentElement == "y":
				#print "    Setting atom y coordinate to %s" % text
				self.currentAtom.y = float(text)
			elif self.currentElement == "z":
				#print "    Setting atom z coordinate to %s" % text
				self.currentAtom.z = float(text)
			elif self.currentElement == "bond":
				#print "    Setting bond text to %s" % text
				self.currentAtom.bonds.append(text)
			elif self.currentElement == "altname":
				#print "    Setting altname text to %s" % text
				self.currentAtom.altname = text
			elif self.currentElement == "dihedral":
				#print "    Setting dihedral text to %s" % text
				self.currentDihedral.atomList = text
			else:
				print "** NOT handling character text:  %s" % text
			

class TopologyResidue:
	""" A class for residue topology information """
	def __init__(self, topology):
		""" Initialize with a Topology object """
		self.name = None
		self.reference = None
		self.titrationStates = []
		self.topology = topology
		self.topology.residues.append(self)
	def __str__(self):
		return self.name
		
class TopologyDihedral:
	""" A class for dihedral topology information.  """
	def __init__(self, parent):
		""" Needs a parent that has a dihedral list. """
		self.parent = parent
		self.parent.dihedrals.append(self)
		self.atomList = None
	def __str__(self):
		return self.atomList
		
class TopologyAtom:
	""" A class for atom topology information """
	def __init__(self, parent):
		""" Needs to be intialized with an upper-level class that contains an atoms array (e.g., TopologyReference
		or TopologyConformerAddition)"""
		self.parent = parent
		self.parent.atoms.append(self)
		self.name = None
		self.x = None
		self.y = None
		self.z = None
		self.bonds = []
		self.altname = None
	def __str__(self):
		return self.name

class TopologyTitrationState:
	""" A class for the titration state of a residue """
	def __init__(self, topologyResidue):
		""" Initialize with a TopologyResidue object """
		self.topologyResidue = topologyResidue
		self.topologyResidue.titrationStates.append(self)
		self.tautomers = []
		self.name = None
	def __str__(self):
		return self.name

class TopologyTautomer:
	""" A class for topology tautomer information """
	def __init__(self, topologyTitrationState):
		""" Initialize with a TopologyTitrationState object """
		self.topologyTitrationState = topologyTitrationState
		self.topologyTitrationState.tautomers.append(self)
		self.conformers = []
		self.name = None
	def __str__(self):
		return self.name		
		
class TopologyConformer:
	""" A class for topology conformer information """
	def __init__(self, topologyTautomer):
		""" Initialize with a TopologyTautomer object """
		self.topologyTautomer = topologyTautomer
		self.topologyTautomer.conformers.append(self)
		self.name = None
		self.conformerAdds = []
		self.conformerRemoves = []
	def __str__(self):
		return self.name
					
class TopologyReference:
	""" A class for the reference structure of a residue """
	def __init__(self, topologyResidue):
		""" Initialize with a TopologyResidue object """
		self.topologyResidue = topologyResidue
		self.topologyResidue.reference = self
		self.atoms = []
		self.dihedrals = []
		
class TopologyConformerAdd:
	""" A class for adding atoms to a conformer """
	def __init__(self, topologyConformer):
		""" Initialize with a TopologyConformer object """
		self.topologyConformer = topologyConformer
		self.topologyConformer.conformerAdds.append(self)
		self.atoms = []
		self.name = None
		self.dihedrals = []

class TopologyConformerRemove:
	""" A class for removing atoms to a conformer """
	def __init__(self, topologyConformer):
		""" Initialize with a TopologyConformer object """
		self.topologyConformer = topologyConformer
		self.topologyConformer.conformerRemoves.append(self)
		self.atoms = []
		self.name = None
		
class Topology:
	""" Contains the structured definitions of residue reference coordinates as well as alternate titration, 
	conformer, and tautomer states.
	"""
	def __init__(self, topologyFile):
		""" Initialize with the topology file reference ready for reading """
		handler = TopologyHandler()
		sax.make_parser()
		sax.parseString(topologyFile.read(), handler)
		self.residues = handler.residues
		
if __name__ == "__main__":
	topologyFile = open(TOPOLOGYPATH, "r")
	topology = Topology(topologyFile)
	
"""
	print "####### SUMMARY ########"
	print "Topology has %d residues:" % len(topology.residues)
	for residue in topology.residues:
		print "-- Residue %s has 1 reference and %d titration states" % (residue.name, len(residue.titrationStates))
		reference = residue.reference
		print "---- Reference %s has %d atoms and %d dihedrals" % (reference, len(reference.atoms), len(reference.dihedrals))
		for atom in reference.atoms:
			print "------ Atom %s" % atom.name
			print "-------- Alt name %s" % atom.altname
			print "-------- Coordinate %g %g %g" % (atom.x, atom.y, atom.z)
		for dihedral in reference.dihedrals:
			print "------ Dihedral %s" % dihedral
		for titrationState in residue.titrationStates:
			print "---- Titration state %s has %d tautomers" % (titrationState.name, len(titrationState.tautomers))
			for tautomer in titrationState.tautomers:
				print "-------- Tautomer %s has %d conformers" % (tautomer.name, len(tautomer.conformers))
				for conformer in tautomer.conformers:
					print "---------- Conformer %s has %d removes" % (conformer.name, len(conformer.conformerRemoves))
					for remove in conformer.conformerRemoves:
						print "------------ Remove %d atoms" % (len(remove.atoms))
						for atom in remove.atoms:
							print "-------------- Atom %s" % (atom.name)
					print "---------- Conformer %s has %d adds" % (conformer.name, len(conformer.conformerAdds))
					for add in conformer.conformerAdds:
						print "------------ Add %d atoms and %d dihedrals" % (len(add.atoms), len(add.dihedrals))
						for atom in add.atoms:
							print "-------------- Atom %s/%s (%g, %g, %g) bonded to %s" % (atom.name, atom.altname, atom.x, atom.y, atom.z, atom.bonds)
						for dihedral in add.dihedrals:
							print "-------------- Dihedral %s" % dihedral
"""
				
		
	
	
