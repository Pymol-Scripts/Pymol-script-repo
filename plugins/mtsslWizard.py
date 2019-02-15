"""
---mtsslWizard: spin_labeling plugin for PyMOL --- 
Author	: Gregor Hagelueken
Date	: July 2015
Version : 1.3
Mail	: hagelueken'at'pc.uni-bonn.de
 
mtsslWizard is a plugin for the PyMOL Molecular Graphics System. 
It allows in silico spin labeling of proteins with  spin labels such as MTSSL in PyMOL. 
Also, distances between ensembles of two spin labels can be calculated and exported.
The program was tested with PyMOL version 1.7.
 
Please cite:
Hagelueken G, Ward R, Naismith JH, Schiemann O. MtsslWizard: In silico Spin-Labeling and Generation of Distance Distributions in PyMOL. 2012. Appl. Mag. Res., accepted for publication.
 
----------------------------------------------------------------------
----------------------------------------------------------------------
 
"""

from __future__ import absolute_import
from __future__ import print_function

import pymol
import string
import threading
import numpy
import scipy.spatial.distance
import random, time, math
import os
from pymol import cmd
from pymol import util
from pymol.wizard import Wizard
from pymol import stored
from operator import itemgetter
import sys

if sys.version_info[0] < 3:
    from Tkinter import Tk
else:
    from tkinter import Tk

default_thoroughness = "thorough search"
default_label = "MTSSL"
default_cutoff = 3.4
default_clashes = 0
default_mode = 'Search'
default_vdwRestraints = "tight"
internalClash_cutoff = 2.5

def __init__(self):
	self.menuBar.addmenuitem('Wizard', 'command',
							 'MtsslWizard',
							 label = 'MtsslWizard',
							 command = lambda s=self : open_wizard())



##########################################################################################
#classes																				 #
##########################################################################################



class MtsslWizard(Wizard):
	def __init__(self):
		#print platform.system()
		Wizard.__init__(self)
		print("")
		print("**************************************************************************************************")
		print("* MtsslWizard by gha.                                                                            *")
		print("* Please remove any solvent or unwanted heteroatoms before using the wizard!                     *")
		print("* You can do this e.g. by issuing 'remove solvent'.                                              *")
		print("**************************************************************************************************")
		print("")
		
		#create array for plots
		try:
			print(stored.plots)
		except:
			stored.plots = []
		
		#create contents of wizard menu
		self.reset()
		self.menu['mode'] = [
									  [1, 'Search','cmd.get_wizard().set_mode("Search")'],
									  [1, 'Measure','cmd.get_wizard().set_mode("Measure")'],
									  [1, 'Distance Map','cmd.get_wizard().set_mode("Distance Map")'],
									  #[1, 'Copy & Move','cmd.get_wizard().set_mode("Copy & Move")']
									  ]
		self.menu['currentLabel'] = [
									  [ 2, '\\559Protein', ''],
									  [1, 'MTSSL','cmd.get_wizard().set_currentLabel("MTSSL")'],
									  [1, 'PROXYL','cmd.get_wizard().set_currentLabel("PROXYL")'],
									  [1, 'DOTA1','cmd.get_wizard().set_currentLabel("DOTA1")'],
									  [1, 'TRITYL','cmd.get_wizard().set_currentLabel("TRITYL")'],
									  #BYSP does not work
									  #[1, 'BYSP','cmd.get_wizard().set_currentLabel("BYSP")'],
									  [1, 'pAcPhe','cmd.get_wizard().set_currentLabel("pAcPhe")'],
									  [1, 'BETA: Rx','cmd.get_wizard().set_currentLabel("Rx")'],
									  [ 2, '\\559DNA', ''],
									  [1, 'URIP','cmd.get_wizard().set_currentLabel("URIP")'],
									  [1, 'C','cmd.get_wizard().set_currentLabel("CLABEL")'],
									  [1, 'TPA','cmd.get_wizard().set_currentLabel("TPA")'],
									  [1, 'Kzn','cmd.get_wizard().set_currentLabel("Kzn")'],
									  ]
# This is the label menu for the DistanceMap Mode
		self.menu['currentLabel1'] = [
									  [ 2, '\\559Protein', ''],
									  [1, 'MTSSL','cmd.get_wizard().set_currentLabel("MTSSL")'],
									  [1, 'PROXYL','cmd.get_wizard().set_currentLabel("PROXYL")'],
									  ]
		self.menu['thoroughness'] = [
									  [1, 'painstaking','cmd.get_wizard().set_thoroughness("painstaking")'],
									  [1, 'thorough search','cmd.get_wizard().set_thoroughness("thorough search")'],
									  [1, 'quick search','cmd.get_wizard().set_thoroughness("quick search")'],
									  ]
		self.menu['vdwRestraints'] = [
									  [ 2, '\\559NOTE:\\559 ', ''],
									  [ 2, '\\559"loose": vdW cutoff 2.5 A, 5 clashes allowed', ''],
									  [ 2, '\\559"tight": vdW cutoff 3.4 A, 0 clashes allowed', ''],
									  [1, 'loose','cmd.get_wizard().set_vdwRestraints("loose")'],
									  [1, 'tight','cmd.get_wizard().set_vdwRestraints("tight")'],
									  ]
		
		self.menu['writeToFile'] = [
									  [1, 'Yes','cmd.get_wizard().set_writeToFile(True)'],
									  [1, 'No','cmd.get_wizard().set_writeToFile(False)']
									  ]

		self.menu['homoOligomerMode'] = [
									  [ 2, '\\955NOTE:\\559Use this to speed up the ', ''],
									  [ 2, '\\559calculation for homooligomers. ', ''],
									  [ 2, '\\559Avoids calculation of redundant maps. ', ''],
									  [1, 'Yes','cmd.get_wizard().set_homoOligomerMode(True)'],
									  [1, 'No','cmd.get_wizard().set_homoOligomerMode(False)']
									  ]

	#some setter and getter functions	  
	def set_vdwRestraints(self, vdwRestraints):
		self.vdwRestraints = vdwRestraints
		if vdwRestraints == "loose":
			cmd.get_wizard().set_cutoff(2.5)
			cmd.get_wizard().set_clashes(5)
		elif vdwRestraints == "tight":
			cmd.get_wizard().set_cutoff(3.4)
			cmd.get_wizard().set_clashes(0)
		self.cmd.refresh_wizard()
	
	def set_rotamers(self,rotamers):
		self.rotamers = rotamers
		self.cmd.refresh_wizard()

	def set_currentLabel(self, currentLabel):
		def mtssl():
			tmp=MtsslLabel("tmp")
			return tmp
		def proxyl():
			tmp=ProxylLabel("tmp")
			return tmp
		def dota1():
			tmp=Dota1Label("tmp")
			return tmp
		def bysp():
			tmp=ByspLabel("tmp")
			return tmp
		def urip():
			tmp=UripLabel("tmp")
			return tmp
		def clabel():
			tmp=CLabel("tmp")
			return tmp
		def pacphe():
			tmp=PAcPheLabel("tmp")
			return tmp
		def tpa():
			tmp=TPA("tmp")
			return tmp
		def trityl():
			tmp=TritylLabel("tmp")
			return tmp
		def rx():
			tmp=RxLabel("tmp")
			return tmp
		def kzn():
			tmp=KznLabel("tmp")
			return tmp
		options = {"MTSSL":mtssl, "PROXYL":proxyl, "DOTA1":dota1, "BYSP":bysp, "URIP":urip, "CLABEL":clabel, "pAcPhe":pacphe, "TPA":tpa, "TRITYL":trityl, "Rx":rx, "Kzn":kzn}
		self.set_vdwRestraints(options[currentLabel]().defaultVdw)
		print(options[currentLabel]().info)
		self.currentLabel = currentLabel
		self.cmd.refresh_wizard()
	
	def set_thoroughness(self,thoroughness):
		self.thoroughness = thoroughness
		self.cmd.refresh_wizard()
		
	def set_writeToFile(self,writeToFile):
		self.writeToFile = writeToFile
		self.cmd.refresh_wizard()
	
	def set_clashGuard(self,clashGuard):
		self.clashGuard = clashGuard
		self.cmd.refresh_wizard()
		
	def set_homoOligomerMode(self,homoOligomerMode):
		self.homoOligomerMode = homoOligomerMode
		self.cmd.refresh_wizard()
		
	def set_cutoff(self,cutoff):
		self.cutoff = cutoff
		self.cmd.refresh_wizard()
		
	def set_clashes(self,clashes):
		self.clashes = clashes
		self.cmd.refresh_wizard()
	
	def set_mode(self, mode):
		self.mode = mode
		self.cmd.refresh_wizard()
			
	def get_prompt(self):
		if self.pick_count == 0 and self.mode == 'Search' and self.currentLabel == 'Rx':
			self.prompt = [ 'Select first anchor point ...']
		if self.pick_count == 1 and self.mode == 'Search' and self.currentLabel == 'Rx':
			self.prompt = [ 'Select second anchor point ...']
		if self.pick_count == 0 and self.mode == 'Search' and self.currentLabel != 'Rx':
			self.prompt = [ 'Select a residue to label...']
		if self.pick_count == 0 and self.mode == 'Measure':
			self.prompt = [ 'Select first label...']
		if self.pick_count == 1 and self.mode == 'Measure':
			self.prompt = [ 'Select second label...']
		if self.running:
			self.prompt = [ 'Running, please wait...' ]
		if self.pick_count == 0 and self.mode == 'Distance Map':
			self.prompt = [ 'Select first object...']
		if self.pick_count == 1 and self.mode == 'Distance Map':
			self.prompt = [ 'Select second object...']
		return self.prompt
	
	def get_panel(self):
		if self.mode == 'Search':
			return [
					[ 1, 'Mtssl Wizard',''],
					[ 3, 'Mode: %s'%self.mode,'mode'],
					[ 3, 'Label: %s'%self.currentLabel,'currentLabel'],
					[ 3, 'Speed: %s'%self.thoroughness,'thoroughness'],
					[ 3, 'vdW restraints: %s'%self.vdwRestraints,'vdwRestraints'],
					[ 2, 'Search conformers!','cmd.get_wizard().run()'],
					[ 2, self.toggleStatesCaption,'cmd.get_wizard().toggle_states()'],
					[ 2, 'Delete last label','cmd.get_wizard().delete_last()'],
					[ 2, 'Reset','cmd.get_wizard().reset()'],
					[ 2, 'Done','cmd.set_wizard()'],
					]
		elif self.mode == 'Measure':
			return [
					[ 1, 'Mtssl Wizard',''],
					[ 3, 'Mode: %s'%self.mode,'mode'],
					[ 3, 'Results to file?: %s'%self.writeToFile,'writeToFile'],
					[ 2, 'Reset','cmd.get_wizard().reset()'],
					[ 2, 'Done','cmd.set_wizard()']
					]
		elif self.mode == 'Distance Map':
			return [
					[ 1, 'MtsslWizard', ''],
					[ 3, 'Mode: %s'%self.mode,'mode'],
					[ 3, 'Label: %s'%self.currentLabel,'currentLabel1'],
					[ 3, 'Homooligomer mode?: %s'%self.homoOligomerMode,'homoOligomerMode'],
					[ 3, 'Results to file?: %s'%self.writeToFile,'writeToFile'],
					[ 2, 'Reset','cmd.get_wizard().reset()'],
					[ 2, 'Done','cmd.set_wizard()']
			]
	#reset wizard to defaults
	def reset(self):
		self.homoOligomerMode = False
		self.running = False
		self.toggleStatesCaption='Toggle states: ON'
		self.start_time=time.time()
		self.conformationList=[]
		self.object_prefix = "mW-"
		self.pick_count = 0
		self.object_count = 0
		self.vdwRestraints="tight"
		self.thoroughness = default_thoroughness
		self.cutoff = default_cutoff
		self.clashes = default_clashes
		self.currentLabel = default_label
		self.residue1_name = None
		self.residue2_name = None
		self.picked_object1 = None
		self.picked_object2 = None
		self.numberOfLabel = 0
		self.mode=default_mode
		self.writeToFile=False
		#self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
		self.label = None
		cmd.set("mouse_selection_mode",1) # set selection mode to residue
		cmd.deselect()
		cmd.unpick()
		cmd.delete(self.object_prefix + "*")
		cmd.delete("sele*")
		cmd.delete("_indicate*")
		cmd.delete("pk*")
		cmd.delete("*_tmp*")
		cmd.refresh_wizard()
 
	def delete_all(self):
		print("Deleting everything...")
		cmd.delete(self.object_prefix + "*")
 
	def delete_last(self):
		try:
			print(self.numberOfLabel)
			if self.numberOfLabel >= 1:
				cmd.delete(self.object_prefix+str(self.numberOfLabel)+"*")
				self.numberOfLabel-=1
		except pymol.CmdException as pmce:
			print(pmce)
	
	def cleanup(self):
		print("Cleaning up...")
		cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
		#self.reset()
		#self.delete_all()
 
	def do_select(self, name):
		try:
			self.do_pick(0)
		except pymol.CmdException as pmce:
			print(pmce)
 
	def do_pick(self, picked_bond):
		if self.mode == 'Search' or self.mode == 'Measure' or self.mode == 'Distance Map':
			#This is to catch repeated selections in Search mode
			if self.pick_count > 0 and self.mode == 'Search' and not self.currentLabel == 'Rx':
				self.pick_count = 0
				cmd.delete("pk*")
				cmd.delete("sele*")
				cmd.delete("*_tmp*")
				cmd.delete("*tmp*")
				cmd.delete("_indicate*")
				cmd.delete("labelEnvironment*")
				self.cmd.refresh_wizard()
				return

			#first click	
			if self.pick_count == 0:	
				self.residue1_name = self.createSelectionMacro("(sele)") 
				# transfer the click selection to a named selection
				cmd.select(self.residue1_name, "(sele)")
				# find the name of the object which contains the selection
				new_name = None
				obj_list = cmd.get_names('objects')
				for object in obj_list:
					if cmd.get_type(object)=="object:molecule":
						if cmd.count_atoms("(%s and %s)"%(object, self.residue1_name)):
							self.picked_object1 = object
							break
				if self.picked_object1 == None:
					print("MtsslWizard: object not found.")
				self.pick_count += 1

				#deselect before next pick
				if self.mode == 'Measure' or self.mode == 'Distance Map' or (self.mode == 'Search' and self.currentLabel == 'Rx'):
					cmd.deselect()
				self.cmd.refresh_wizard()

			#second click
# 			if self.mode == 'Search' and self.currentLabel == 'Rx':
# 				print "Test!"
# 				print self.pick_count
# 			if self.mode == 'Measure' or self.mode == 'Distance Map' or (self.mode == 'Search' and self.currentLabel == 'Rx'):
# 				print "Another Test"
			elif self.pick_count == 1 and (self.mode == 'Measure' or self.mode == 'Distance Map' or (self.mode == 'Search' and self.currentLabel == 'Rx')):
				self.residue2_name = self.createSelectionMacro("(sele)")
# 				print self.residue2_name
# 				print "Hallo, Test!"
				# transfer the click selection to a named selection
				cmd.select(self.residue2_name, "(sele)")
				# find the name of the object which contains the selection
				new_name = None
				obj_list = cmd.get_names('objects')
				for object in obj_list:
					if cmd.get_type(object)=="object:molecule":
						if cmd.count_atoms("(%s and %s)"%(object, self.residue2_name)):
							self.picked_object2 = object
							break

				if self.picked_object2 == None:
					print("MtsslWizard: object not found.")
				self.pick_count += 1
				#deselect before next pick
				if self.mode == "Distance Map" or self.mode == "Measure" or (self.mode == "Search" and self.currentLabel == "Rx"):
					self.run()
				cmd.deselect()
				#print self.picked_object1
				#print self.picked_object2
		
			if self.mode == 'Search':
				self.numberOfLabel += 1
	
	#def delay_launch(self):
	#	self.running = True
	#	self.cmd.refresh_wizard()
	#	cmd.draw()
	#	self.cmd.feedback("disable","all","everything")
	#	self.cmd.feedback("enable","python","output")
	#	t = threading.Thread(target=self.run)
	#	t.setDaemon(1)
	#	t.start()
	#	#time.sleep(10)
	#	#self.run()
	#	self.running = False
	#	self.cmd.refresh_wizard()


	def setLabel(self, name):
		if self.currentLabel == "MTSSL":
			self.label = MtsslLabel(name)
			
		elif self.currentLabel == "PROXYL":
			self.label = ProxylLabel(name)
			
		elif self.currentLabel == "DOTA1":
			self.label = Dota1Label(name)
			
		elif self.currentLabel == "URIP":
			self.label = UripLabel(name)
			
		elif self.currentLabel == "CLABEL":
			self.label = CLabel(name)
			
		elif self.currentLabel == "BYSP":
			self.label = ByspLabel(name)
			
		elif self.currentLabel == "pAcPhe":
			self.label = PAcPheLabel(name)
			
		elif self.currentLabel == "TPA":
			self.label = TPA(name)
		
		elif self.currentLabel == "TRITYL":
			self.label = TritylLabel(name)
		
		elif self.currentLabel == "Rx":
			self.label = RxLabel(name)
		
		elif self.currentLabel == "Kzn":
			self.label = KznLabel(name)
			
		
##########################################################################################
#Run																					 #
##########################################################################################
	def run(self):
		self.conformationList = []
		my_view = cmd.get_view()
##########################################################################################
#Search																					 #
##########################################################################################
		
		if self.mode == 'Search':
			#generate the label and superpose onto selected position
			if self.currentLabel == "MTSSL":
				self.label = MtsslLabel("mtssl_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "PROXYL":
				self.label = ProxylLabel("proxyl_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "DOTA1":
				self.label = Dota1Label("dota1_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "URIP":
				self.label = UripLabel("urip_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "CLABEL":
				self.label = CLabel("C_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "BYSP":
				self.label = ByspLabel("bysp_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "pAcPhe":
				self.label = PAcPheLabel("pAcPhe_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "TPA":
				self.label = TPA("tpa_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "TRITYL":
				self.label = TritylLabel("TRITYL_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "Rx":
				self.label = RxLabel("Rx"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "Kzn":
				self.label = KznLabel("Kzn"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			print(self.residue1_name, self.label.identifier)
		
		if self.pick_count == 1 and self.mode == 'Search':
			print("\n\n\nNew run:\n")	
			#attach identifier for identification in pymol
			#cmd.set_name(self.residue1_name, "%s_%s" %(self.residue1_name, self.label.identifier))
			#self.residue1_name = "%s_%s" %(self.residue1_name, self.label.identifier)
			
			print("Attempting superposition...")
			if not self.superpose():
				print("Superposition does not work.")
				print("Possible reasons:")
				print("1) Glycine? Mutate to Ala first.")
				print("2) Trying to attach DNA label to Protein or vice versa?")
				if len(self.label.errorMessage) > 0:
					print("3) %s" %self.label.errorMessage)
				self.cleanupAfterRun(my_view)
				return
			else:
				print("Superposition worked!")
			
			#prepare movingAtoms array of label, put into correct order...
			stored.movingAtoms = []
			for i in range (0, len(self.label.atomNames)):
				xyz = cmd.get_model("%s & name %s" %(self.label.pymolName, self.label.atomNames[i] ), 1).get_coord_list()
				stored.movingAtoms.extend(xyz)
			self.label.movingAtoms=numpy.array(stored.movingAtoms)
			
			#create object with only the atoms around the label to speed everything up 
			#cmd.color ("red", "%s &! %s within %f of %s" %(self.picked_object1, self.residue1_name, self.label.radius, self.label.pymolName))
			protein ="%s &! %s within %f of %s" %(self.picked_object1, self.residue1_name, self.label.radius, self.label.pymolName)
			cmd.create ("labelEnvironment", "byres "+protein)
			stored.environmentAtoms = []
			cmd.iterate_state(1, protein, 'stored.environmentAtoms.append((x,y,z))')
			environmentAtoms = numpy.array(stored.environmentAtoms)
	
			result = self.fastMtsslify(environmentAtoms)
			#only switch on snuggly fit search for "painstaking"
			if self.thoroughness == "painstaking":
				self.createSnugglyFitConformations()		
			print("")
			print("Found: %i in %i tries." %(result[0], result[1]))
			if result[0] > 0 and result[0] <= 10 and self.vdwRestraints == "tight" and not self.currentLabel == "CLABEL":
				print("The number of conformations is very small!")
				print("Consider switching vdW restraints to 'loose' for this position!")
			print("Done!")
			self.finalCosmetics(result[0])
			self.cleanupAfterRun(my_view)
		
		#Rx label
		if self.pick_count == 2 and self.mode == 'Search' and self.currentLabel == 'Rx':
			print("\n\n\nNew run:\n")
			ca1 = numpy.array(cmd.get_model(self.residue1_name + " & name CA", 1).get_coord_list()[0])
			ca2 = numpy.array(cmd.get_model(self.residue2_name + " & name CA", 1).get_coord_list()[0])
			try:
				cb1 = numpy.array(cmd.get_model(self.residue1_name + " & name CB", 1).get_coord_list()[0])
			except:
				cb1 = ca1
			try:
				cb2 = numpy.array(cmd.get_model(self.residue2_name + " & name CB", 1).get_coord_list()[0])
			except:
				cb2 = ca2
			environmentatoms = numpy.array(cmd.get_model("%s within 10 of %s or %s and not (%s or %s)" %(self.picked_object1, self.residue1_name, self.residue2_name, self.residue1_name, self.residue2_name), 1).get_coord_list())
			averageConeCoordinate1, numberOfConeAtoms1, solutions1 = self.calculateCone(ca1, cb1, environmentatoms, numberOfAtoms=4000)
			averageConeCoordinate2, numberOfConeAtoms2, solutions2 = self.calculateCone(ca2, cb2, environmentatoms, numberOfAtoms=4000)
			distances1 = quick_map(solutions1, cb2)
			distances2 = quick_map(solutions2, cb1)
			#print distances
			indices1 = numpy.where(numpy.any(distances1 > 8, axis=1))
			#print indices1
			indices2 = numpy.where(numpy.any(distances2 > 8, axis=1))
			#print indices2
			solutions1 = numpy.delete(solutions1, indices1, 0)
			solutions2 = numpy.delete(solutions2, indices2, 0)
			solutions = numpy.concatenate((solutions1, solutions2))
			#self.createRxRotamer(solutions)
			for idx, solution in enumerate(solutions):
				#print solution
				self.label.movingAtoms = solution
				#print self.label.movingAtoms
				self.createRotamer(idx)
			self.cleanupAfterRun(my_view)
			if len(solutions) == 0:
				print("Did not find any possible Rx positions. Are the two residues too far apart?")
			else:
				#print "Found %i possible Rx positions." %len(solutions)
				self.finalCosmetics(len(solutions))
			


##########################################################################################
#Distance map																			 #
##########################################################################################
		
		elif self.pick_count == 2 and self.mode == 'Distance Map':
			self.setLabel("dummy")
			
			#if both picked objects are identical
			if self.picked_object1 == self.picked_object2:
				labelPositions = []
				stored.residueList = []
				chains = cmd.get_chains(self.picked_object1)
				if self.homoOligomerMode:
					chains_1 = chains[0]
				else:
					chains_1 = chains
				#iterate over chains
				for chain_1 in chains_1:
					for chain_2 in chains:
						chain_1_residues = []
						cmd.iterate("%s & chain %s & name CA" %(self.picked_object1, chain_1), "stored.residueList.append((model, chain, resi,resn))")
						chain_1_residues = stored.residueList
						stored.residueList = []
						chain_2_residues = []
						cmd.iterate("%s & chain %s & name CA" %(self.picked_object1, chain_2), "stored.residueList.append((model, chain, resi,resn))")
						chain_2_residues = stored.residueList
						stored.residueList = []
						#print stored.residueList
						#return
						chain_1_labelPositions = []
						chain_2_labelPositions = []
						chain_1_cBetaPositions = []
						chain_2_cBetaPositions = []
						chain_1_numberOfConeAtoms = []
						chain_2_numberOfConeAtoms = []
						
						#iterate over residues in chain 1 and calculate approximate label position
						for residue in chain_1_residues:
							selectionString = "%s & chain %s & resi %s" %(residue[0], residue[1], residue[2])
							#print selectionString
							ca = numpy.array(cmd.get_model(selectionString + " & name CA", 1).get_coord_list()[0])
							try:
								cb = numpy.array(cmd.get_model(selectionString + " & name CB", 1).get_coord_list()[0])
							except:
								#print "Glycine: %s" %selectionString
								#print "\t\tUsing C-alpha as label position..."
								cb = ca
							#try:
							environmentatoms = numpy.array(cmd.get_model("(%s within 10 of %s) &! (%s)" %(self.picked_object1, selectionString, selectionString), 1).get_coord_list())
							averageConeCoordinate, numberOfConeAtoms, solutions = self.calculateCone(ca, cb, environmentatoms)
# 							except Exception, e:
# 								print "Could not label: %s" %selectionString
# 								print "\t\tNot on surface?"
# 								print "\t\tUsing C-alpha as label position..."
# 								print "Error message: %s" %e
# 								averageConeCoordinate = ca # + numpy.random.rand(3)
							chain_1_labelPositions.append(averageConeCoordinate)
							chain_1_cBetaPositions.append(cb)
							chain_1_numberOfConeAtoms.append(numberOfConeAtoms)
						
						#iterate over residues in chain 2 and calculate approximate label position
						for residue in chain_2_residues:
							selectionString = "%s & chain %s & resi %s" %(residue[0], residue[1], residue[2])
							#print selectionString
							ca = numpy.array(cmd.get_model(selectionString + " & name CA", 1).get_coord_list()[0])
							try:
								cb = numpy.array(cmd.get_model(selectionString + " & name CB", 1).get_coord_list()[0])
							except:
								#print "Glycine: %s" %selectionString
								#print "		Using C-alpha as label position..."
								cb = ca
							try:
								environmentatoms = numpy.array(cmd.get_model("(%s within 10 of %s) &! (%s)" %(self.picked_object2, selectionString, selectionString), 1).get_coord_list())
								averageConeCoordinate, numberOfConeAtoms, solutions = self.calculateCone(ca, cb, environmentatoms)
							except:
								#print "Could not label: %s" %selectionString
								#print "		Not on surface?"
								#print "		Using C-alpha as label position..."
								averageConeCoordinate = ca # + numpy.random.rand(3)
							chain_2_labelPositions.append(averageConeCoordinate)
							chain_2_cBetaPositions.append(cb)
							chain_2_numberOfConeAtoms.append(numberOfConeAtoms)
						
						#calculate the distance matrix and its diagonal
						distances = numpy.nan_to_num(quick_map(chain_1_labelPositions, chain_2_labelPositions))
						diagonal = numpy.diagonal(distances)
						
						#write out the data and plots
						fileName = "%s_%s-%s_distanceMatrix" %(self.picked_object1, chain_1, chain_2)
						if self.writeToFile:
							numpy.savetxt(fileName+".txt", distances)
						offset1 = int(chain_1_residues[0][2])
						offset2 = int(chain_2_residues[0][2])
						x = numpy.linspace(offset1, distances.shape[0]+offset1, distances.shape[0]+1)
						y = numpy.linspace(offset2, distances.shape[1]+offset2, distances.shape[1]+1)
						z = numpy.transpose(distances)
						#plot=Plotter()
						#plot.writeToFile = self.writeToFile
						#try:
						#	plot.plotDistanceMap(x, y, z, fileName+".png")
							
						#except:
						#	print "Could not plot map."
						#for mtsslPlot
						plotDictionary = self.makeGraphDataDictionary(fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, z)
						stored.plots.append(plotDictionary)
						
						#plot.plotOffsetMap(x, y, numpy.transpose(distances), fileName+".png")
						if self.writeToFile:
							numpy.savetxt("%s_%s-%s_distanceMatrix.txt" %(self.picked_object1,chain_1,chain_2), distances)
						
						#Plot accessibility for both chains
						#chain 1
						fileName = "%s_%s_accessibilty" %(self.picked_object1, chain_1)
						x = numpy.linspace(offset1, distances.shape[0]+offset1-1, distances.shape[0])
						y = numpy.array(chain_1_numberOfConeAtoms)
						#plot.plotAccessibility(x, y, fileName+".png")
						if self.writeToFile:
							numpy.savetxt("%s_%s_accessibility.txt" %(self.picked_object1,chain_1), numpy.column_stack((x,numpy.array(chain_1_numberOfConeAtoms))))
						#for mtsslPlot
						plotDictionary = self.makeGraphDataDictionary(fileName, "AccessibilityPlot", "Number of Residue", "Relative Accessibility", x, y, 0)
						stored.plots.append(plotDictionary)
						#chain 2
						fileName = "%s_%s_accessibilty" %(self.picked_object1, chain_2)
						x = numpy.linspace(offset2, distances.shape[1]+offset2-1, distances.shape[1])
						y = numpy.array(chain_2_numberOfConeAtoms)
						#plot.plotAccessibility(x, y, fileName+".png")
						if self.writeToFile:
							numpy.savetxt("%s_%s_accessibility.txt" %(self.picked_object1,chain_2), numpy.column_stack((x,numpy.array(chain_2_numberOfConeAtoms))))
						#for mtsslPlot
						plotDictionary = self.makeGraphDataDictionary(fileName, "AccessibilityPlot", "Number of Residue", "Relative Accessibility", x, y, 0)
						stored.plots.append(plotDictionary)
						
						
						#write out diagonal only if both chains have the same number of residues
						if distances.shape[0] == distances.shape[1]:
							fileName = "%s_%s-%s_diagonal" %(self.picked_object1,chain_1,chain_2)
							if self.writeToFile:
								numpy.savetxt("%s_%s-%s_diagonal.txt" %(self.picked_object1,chain_1,chain_2), diagonal)
							offset1 = int(chain_1_residues[0][2])
							#10 110 100
							offset2 = int(chain_2_residues[0][2])
							x = numpy.linspace(offset1, distances.shape[0]+offset1-1, distances.shape[0])
							y = diagonal
							#plot=Plotter()
							if self.writeToFile:
								numpy.savetxt(fileName+".txt", numpy.column_stack((x, diagonal)))
							try:
								#plot.plotDistances(x, y, fileName+".png")
								#for mtsslPlot
								plotDictionary = self.makeGraphDataDictionary(fileName, "DistancePlot", "Number of Residue", "Distance (Angstrom)", x, y, 0)
								stored.plots.append(plotDictionary)
							except:
								print("Could not plot map!")
						print("Still working...")
			
			#if different objects were picked 
			elif self.picked_object1 != self.picked_object2:
				labelPositions_object1 = []
				labelPositions_object2 = []
				cbetaPositions_object1 = []
				cbetaPositions_object2 = []
				numberOfConeAtoms_object1 = []
				numberOfConeAtoms_object2 = []
				
				stored.residueList = []
				chains_object1 = cmd.get_chains(self.picked_object1)
				chains_object2 = cmd.get_chains(self.picked_object2)
				#print self.picked_object1, chains_object1
				#print self.picked_object2, chains_object2
				#iterate over chains in object 1
				for chain_1 in chains_object1:
					labelPositions_object1 = []
					cbetaPositions_object1 = []
					chain_1_residues = []
					stored.residueList = []
					cmd.iterate("%s & chain %s & name CA" %(self.picked_object1, chain_1), "stored.residueList.append((model, chain, resi,resn))")
					chain_1_residues = stored.residueList
					stored.residueList = []
					#iterate over residues in chain 1 and calculate approximate spin label position
					for residue in chain_1_residues:
						selectionString = "%s & chain %s & resi %s" %(residue[0], residue[1], residue[2])
						#print selectionString
						ca = numpy.array(cmd.get_model(selectionString + " & name CA", 1).get_coord_list()[0])
						try:
							cb = numpy.array(cmd.get_model(selectionString + " & name CB", 1).get_coord_list()[0])
						except:
							#print "Glycine: %s" %selectionString
							#print "\t\tUsing C-alpha as label position..."
							cb = ca
						try:
							environmentatoms = numpy.array(cmd.get_model("(%s within 10 of %s) &! (%s)" %(self.picked_object1, selectionString, selectionString), 1).get_coord_list())
							averageConeCoordinate, numberOfConeAtoms, solutions = self.calculateCone(ca, cb, environmentatoms)
						except:
							#print "Could not label: %s" %selectionString
							#print "\t\tNot on surface?"
							#print "\t\tUsing C-alpha as label position..."
							averageConeCoordinate = ca # + numpy.random.rand(3)
						labelPositions_object1.append(averageConeCoordinate)
						cbetaPositions_object1.append(cb)
						numberOfConeAtoms_object1.append(numberOfConeAtoms)
					if cmd.count_atoms(self.picked_object1) == 1:
						labelPositions_object1.append(numpy.array(cmd.get_model(self.picked_object1, 1).get_coord_list()[0]))
						
					#calculate intra chain distances for current chain in object 1
					intraDistancesObject1 = quick_map(labelPositions_object1, labelPositions_object1)
					cbetaDistancesObject1 = quick_map(cbetaPositions_object1, cbetaPositions_object1)
					
					#iterate over chains in object 2
					for chain_2 in chains_object2:
						labelPositions_object2 = []
						cbetaPositions_object2 = []
						chain_2_residues = []
						stored.residueList = []
						cmd.iterate("%s & chain %s & name CA" %(self.picked_object2, chain_2), "stored.residueList.append((model, chain, resi,resn))")
						chain_2_residues = stored.residueList
						stored.residueList = []
						for residue in chain_2_residues:
							selectionString = "%s & chain %s & resi %s" %(residue[0], residue[1], residue[2])
							#print selectionString
							ca = numpy.array(cmd.get_model(selectionString + " & name CA", 1).get_coord_list()[0])
							try:
								cb = numpy.array(cmd.get_model(selectionString + " & name CB", 1).get_coord_list()[0])
							except:
								#print "Glycine: %s" %selectionString
								#print "		Using C-alpha as label position..."
								cb = ca
							try:
								environmentatoms = numpy.array(cmd.get_model("(%s within 10 of %s) &! (%s)" %(self.picked_object2, selectionString, selectionString), 1).get_coord_list())
								averageConeCoordinate, numberOfConeAtoms, solutions = self.calculateCone(ca, cb, environmentatoms)
							except:
								#print "Could not label: %s" %selectionString
								#print "		Not on surface?"
								#print "		Using C-alpha as label position..."
								averageConeCoordinate = ca # + numpy.random.rand(3)
							labelPositions_object2.append(averageConeCoordinate)
							cbetaPositions_object2.append(cb)
							numberOfConeAtoms_object2.append(numberOfConeAtoms)
						if cmd.count_atoms(self.picked_object2) == 1:
							labelPositions_object2.append(numpy.array(cmd.get_model(self.picked_object2, 1).get_coord_list()[0]))
						#calculate intra chain distances for for current chain in object 2
						intraDistancesObject2 = quick_map(labelPositions_object2, labelPositions_object2)
						cbetaDistancesObject2 = quick_map(cbetaPositions_object2, cbetaPositions_object2)
						calculatedDifferenceDistanceMatrix = False
						#calculate inter chain distances for object 1-2
						interDistancesObject1_2 = quick_map(labelPositions_object1, labelPositions_object2)
						interDistancesObject2_1 = quick_map(labelPositions_object2, labelPositions_object1)
						
						#calculate difference distance matrix
						try:
							differences = numpy.abs(intraDistancesObject1-intraDistancesObject2)
							cbetaDifferences = numpy.abs(cbetaDistancesObject1-cbetaDistancesObject2)
							calculatedDifferenceDistanceMatrix = True
						except:
							print(len(intraDistancesObject1))
							print(len(intraDistancesObject2))
							print(len(cbetaDistancesObject1))
							print(len(cbetaDistancesObject2))
							print("\nCannot subtract matrices. Please check if both objects have the same number of residues!\nRemove alternate conformations.\n")
							print("Try:")
							print("remove not (alt ''+A)")
							print("alter all, alt=''")
						
						#plot and write distance maps
						if cmd.count_atoms(self.picked_object1) > 1 and cmd.count_atoms(self.picked_object2) > 1:
							
							#intermolecular map 1
							fileName = "%s_%s_%s-%s_distanceMatrix" %(self.picked_object1, self.picked_object2, chain_1, chain_2)
							if self.writeToFile:
								numpy.savetxt(fileName+".txt", interDistancesObject1_2)
							try:
								offset1 = int(chain_1_residues[0][2])
							except:
								offset1 = 0
							try:
								offset2 = int(chain_2_residues[0][2])
							except:
								offset2 = 0
							x = numpy.linspace(offset1, interDistancesObject1_2.shape[0]+offset1, interDistancesObject1_2.shape[0]+1)
							y = numpy.linspace(offset2, interDistancesObject1_2.shape[1]+offset2, interDistancesObject1_2.shape[1]+1)
							z = numpy.transpose(interDistancesObject1_2)
							#plot=Plotter()
							#plot.writeToFile = self.writeToFile
							try:
								#plot.plotDistanceMap(x, y, z, fileName+".png")
								plotDictionary = self.makeGraphDataDictionary(fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, z)
								stored.plots.append(plotDictionary)
							except:
								print("Could not plot map!")
						
							#intramolecular map 1
							fileName = "%s_%s_distanceMatrix.txt" %(self.picked_object1,  chain_1)
							if self.writeToFile:
								numpy.savetxt(fileName+".txt", intraDistancesObject1)
							try:
								offset1 = int(chain_1_residues[0][2])
							except:
								offset1 = 0
							try:
								offset2 = int(chain_2_residues[0][2])
							except:
								offset2 = 0
							x = numpy.linspace(offset1, intraDistancesObject1.shape[0]+offset1, intraDistancesObject1.shape[0]+1)
							y = numpy.linspace(offset2, intraDistancesObject1.shape[1]+offset2, intraDistancesObject1.shape[1]+1)
							z = numpy.transpose(intraDistancesObject1)
							try:
								#plot.plotDistanceMap(x, y, z, fileName+".png")
								plotDictionary = self.makeGraphDataDictionary(fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, z)
								stored.plots.append(plotDictionary)
							except:
								print("Could not plot map!")
						
							#intramolecular map 2
							fileName = "%s_%s_distanceMatrix.txt" %(self.picked_object2,  chain_2)
							if self.writeToFile:
								numpy.savetxt(fileName+".txt", intraDistancesObject2)
							try:
								offset1 = int(chain_1_residues[0][2])
							except:
								offset1 = 0
							try:
								offset2 = int(chain_2_residues[0][2])
							except:
								offset2 = 0
							x = numpy.linspace(offset1, intraDistancesObject2.shape[0]+offset1, intraDistancesObject2.shape[0]+1)
							y = numpy.linspace(offset2, intraDistancesObject2.shape[1]+offset2, intraDistancesObject2.shape[1]+1)
							z = numpy.transpose(intraDistancesObject2)
							
							try:
								#plot.plotDistanceMap(x, y, z, fileName+".png")
								plotDictionary = self.makeGraphDataDictionary(fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, z)
								stored.plots.append(plotDictionary)
							except:
								print("Could not plot map!")
							
							
							#difference distance map
							if calculatedDifferenceDistanceMatrix:
								fileName = "%s_%s_%s-%s_differenceDistanceMatrix" %(self.picked_object1, self.picked_object2, chain_1, chain_2)
								if self.writeToFile:
									numpy.savetxt(fileName+".txt", differences)
								try:
									offset1 = int(chain_1_residues[0][2])
								except:
									offset1 = 0
								try:
									offset2 = int(chain_2_residues[0][2])
								except:
									offset2 = 0
								#diagonal = numpy.diagonal(differences)
								#print diagonal
								x = numpy.linspace(offset1, differences.shape[0]+offset1, differences.shape[0]+1)
								y = numpy.linspace(offset2, differences.shape[1]+offset2, differences.shape[1]+1)
								
								z1 = numpy.transpose(differences)
								z2 = numpy.transpose(cbetaDifferences)
								try:
									#plot.plotOffsetMap(x, y, numpy.transpose(differences), fileName+".png")
									plotDictionary = self.makeGraphDataDictionary(fileName, "DifferenceMap", "Number of Residue", "Number of Residue", x, y, z1)
									stored.plots.append(plotDictionary)
									#plot.plotOffsetMap(x, y, numpy.transpose(cbetaDifferences), fileName+"-cbeta.png")
									plotDictionary = self.makeGraphDataDictionary(fileName+"-cbeta", "DifferenceMap", "Number of Residue", "Number of Residue", x, y, z2)
									stored.plots.append(plotDictionary)
									#plotDictionary = self.makeGraphDataDictionary(fileName+"_diagonal", "DistancePlot", "Number of Residue", "Distance (Angstrom)", x, diagonal, 0)
									#stored.plots.append(plotDictionary)
								except:
									print("Could not plot map!")
						
						#if one object was a single atom, make a scatter plot
						if cmd.count_atoms(self.picked_object1) == 1 or cmd.count_atoms(self.picked_object2) == 1:
							
							try:
								offset1 = int(chain_1_residues[0][2])
							except:
								offset1 = 0
							try:
								offset2 = int(chain_2_residues[0][2])
							except:
								offset2 = 0
							
							if cmd.count_atoms(self.picked_object1) == 1:
								x = numpy.linspace(offset2, interDistancesObject2_1.shape[0]+offset2-1, interDistancesObject2_1.shape[0])
								y = interDistancesObject2_1
							if cmd.count_atoms(self.picked_object2) == 1:
								x = numpy.linspace(offset1, interDistancesObject1_2.shape[0]+offset1-1, interDistancesObject1_2.shape[0])
								y = interDistancesObject1_2
							
							fileName = "%s_%s_%s-%s_distancePlot" %(self.picked_object1, self.picked_object2, chain_1, chain_2)
							if self.writeToFile:
								numpy.savetxt(fileName+".txt", numpy.column_stack((x, y)))
							try:
								#plot.plotDistances(x, y, fileName+".png")
								plotDictionary = self.makeGraphDataDictionary(fileName, "DistancePlot", "Number of Residue", "Distance (Angstrom)", x, y, 0)
								stored.plots.append(plotDictionary)
							except:
								print("Could not plot map!")
						
			self.cleanupAfterRun(my_view)
			print("Done!")
							

##########################################################################################
#Measure																				 #
##########################################################################################
		elif self.pick_count == 2 and self.mode == 'Measure':
			print("\n\n\nDistance calculation:\n")
			print("The dashed lines are the c-beta distance (green),\nand the distance between the geometric averages\nof the two ensembles (yellow).\n")
			print("The following statistics refer to the distribution\nof the individual distances between all conformers (may take a while):\n")
			#find out what the selections are
			stored.label1 = []
			stored.label2 = []
			stored.atomNames1 = []
			stored.atomNames2 = []
			cmd.iterate(self.residue1_name, 'stored.atomNames1.append(name)')
			cmd.iterate(self.residue2_name, 'stored.atomNames2.append(name)')
			
			#This is to find the spin location of the selected labels
			def mtssl():
				tmp=MtsslLabel("tmp")
				return tmp.spinLocation
			def proxyl():
				tmp=ProxylLabel("tmp")
				return tmp.spinLocation
			def dota1():
				tmp=Dota1Label("tmp")
				return tmp.spinLocation
			def bysp():
				tmp=ByspLabel("tmp")
				return tmp.spinLocation
			def urip():
				tmp=UripLabel("tmp")
				return tmp.spinLocation
			def clabel():
				tmp=CLabel("tmp")
				return tmp.spinLocation
			def pacphe():
				tmp=PAcPheLabel("tmp")
				return tmp.spinLocation
			def tpa():
				tmp=TPA("tmp")
				return tmp.spinLocation
			def trityl():
				tmp=TritylLabel("tmp")
				return tmp.spinLocation
			def rx():
				tmp=RxLabel("tmp")
				return tmp.spinLocation
			def kzn():
				tmp=KznLabel("tmp")
				return tmp.spinLocation
			options = {"M-T-S-S-L":mtssl, "P-R-O-X-Y-L":proxyl, "D-O-T-A-1":dota1, "B-Y-S-P":bysp, "U-R-I-P":urip, "C-L-A-B-E-L":clabel, "p-A-c-P-h-e":pacphe, "T-P-A":tpa, "T-R-I-T-Y-L":trityl, "R-X":rx, "K-Z-N":kzn}
			
			#Decide if only the spin location (for labels) or all atoms of the selection are used
			if self.picked_object1.split('_')[-1] in options:
				spinLocation = options[self.picked_object1.split('_')[-1]]()
				cmd.iterate_state(0, "%s & name %s" %(self.residue1_name, spinLocation), 'stored.label1.append((x,y,z))')
			else:	
				cmd.iterate_state(0, self.residue1_name, 'stored.label1.append((x,y,z))')
			if self.picked_object2.split('_')[-1] in options:
				spinLocation = options[self.picked_object2.split('_')[-1]]()
				cmd.iterate_state(0, "%s & name %s" %(self.residue2_name, spinLocation), 'stored.label2.append((x,y,z))')
			else:
				cmd.iterate_state(0, self.residue2_name, 'stored.label2.append((x,y,z))')
			
			atoms1=numpy.array(stored.label1)
			atoms2=numpy.array(stored.label2)
			#Calculate the distances
			dist=quick_dist2(atoms1, atoms2)
			
			#create pseudoatom at average coordinate of each ensemble and display the distance between them
			avgAtoms1=numpy.average(atoms1,axis=0)
			avgAtoms2=numpy.average(atoms2,axis=0)
			self.createPseudoatom (avgAtoms1, "tmp_average1", 1)
			self.createPseudoatom (avgAtoms2, "tmp_average2", 1)
			cmd.distance(self.object_prefix+"avg","tmp_average1 & name PS1","tmp_average2 & name PS1")
			cmd.delete("tmp_average1")
			cmd.delete("tmp_average2")
			
			#cbeta distance if cbeta is present in both selections
			cBeta = []
			if any("CB" in atom for atom in stored.atomNames1) and any("CB" in atom for atom in stored.atomNames2):
				cmd.distance(self.object_prefix+"cBeta", self.residue1_name+" & name CB",self.residue2_name+" & name CB")
				cmd.set("dash_color", "green", self.object_prefix+"cBeta")
				stored.label1 = []
				stored.label2 = []
				cmd.iterate_state(1, self.residue1_name+" & name CB", 'stored.label1.append((x,y,z))')
				cmd.iterate_state(1, self.residue2_name+" & name CB", 'stored.label2.append((x,y,z))')
				atoms1=numpy.array(stored.label1)
				atoms2=numpy.array(stored.label2)
				cBeta=quick_dist2(atoms1, atoms2)
			
			#output of distances and histogram to clipboard or file
			numpy.set_printoptions(threshold=10000000, precision = 2, suppress = True)
			#copy to clipboard
			#create envelope plot
			distListForOutput = numpy.copy(dist)
			if len(dist) < 100:
				#fill up with 999999s which can be easily removed from output
				distListForOutput.resize(100)
				distListForOutput[len(dist):]=9999999
			histogram=numpy.histogram(distListForOutput, numpy.arange(100))
			envelopePlot = numpy.zeros((len(distListForOutput),2))
			envelopePlot[0:99] = numpy.column_stack((histogram[1][0:len(histogram[1])-1], histogram[0]))
			#put point in mid of bin
			envelopePlot[:,0] += 0.5 
			normEnvelopePlot = numpy.copy(envelopePlot)
			normEnvelopePlot[:,1] = normEnvelopePlot[:,1]/numpy.amax(histogram[0])
			#combine dist and histogram to single array before output
			output=numpy.column_stack((distListForOutput, envelopePlot, normEnvelopePlot[:,1]))
			outputStr=numpy.array_str(output)
			outputStr=outputStr.replace("[", "")
			outputStr=outputStr.replace("]", "")
			outputStr=outputStr.replace("9999999","")
			graphtitle = "%s-%s" %(self.residue1_name, self.residue2_name)
			plotDictionary = self.makeGraphDataDictionary (graphtitle, "DistanceDistribution", "Distance (Angstrom)", "Relative Probability", output[:,1], output[:,3], 0)
			stored.plots.append(plotDictionary)
			print("Distribution plot added to memory. Check it with mtsslPlotter.")
			
			#Copy to clipboard
			self.copyStringToClipboard(outputStr)

			#Write to file
			if self.writeToFile:
				try:
					filename = "%s-%s" %(self.residue1_name, self.residue2_name)
					numpy.savetxt(filename, output, delimiter='\t')
					print("Written to file:")
					print("%s/%s" %(os.getcwd(), filename))
				except:
					print("Writing to file failed!")
					
			print(calculateStatistics2(dist))
			if len(cBeta) > 0:
				print("Cbeta distance: %3.1f" %cBeta[0])
			self.cleanupAfterRun(my_view)
			print("Done.")
	
##########################################################################################
#various methods																				 #
##########################################################################################
	
	def makeGraphDataDictionary (self, title, type, xTitle, yTitle, xData, yData, zData):
		graphData = {}
		graphData['Title'] = title
		graphData['Type'] = type
		graphData['xTitle'] = xTitle
		graphData['yTitle'] = yTitle
		graphData['xData'] = xData
		graphData['yData'] = yData
		graphData['zData'] = zData
		return graphData
	
	def calculateApproxLabelPosition(self, ca, cb, n, c, name):
		#print ca, cb
		caCbVector = (cb-ca)/numpy.linalg.norm(cb-ca)
		#print caCbVector
		#createPseudoatom(caCbVector, "caCb", 1)
		labelVector = caCbVector * 8.0
		labelPosition = ca + labelVector
		#self.createPseudoatom(labelPosition, name, 1)
		return labelPosition
	
	def copyStringToClipboard(self, string):
		#try:
		#	r = pymol._ext_gui.root
		#	r.clipboard_clear()
		#	r.clipboard_append(string)
		#	print "Copied to clipboard."
		#	return
		#except:
		#	pass
		try:
			import pyperclip
			pyperclip.copy(string)
			print("Copied to clipboard.")
			return
		except:
			pass
		try:
			import xerox
			xerox.copy(string)
			print("Copied to clipboard.")
			return
		except:
			pass
		print("Copy to clipboard failed. Try to install either the 'pyperclip' or 'xerox' module for Python.")
	
	def cleanupAfterRun(self, my_view):
		self.pick_count = 0
		cmd.delete("pk*")
		cmd.delete("sele*")
		cmd.delete("*_tmp*")
		cmd.delete("*tmp*")
		cmd.delete("_indicate*")
		cmd.delete("labelEnvironment*")
		self.cmd.refresh_wizard()
		cmd.set_view(my_view)
	
	def toggle_states(self):
		if cmd.get("all_states")=='on':
			self.toggleStatesCaption='Toggle states: OFF'
			cmd.set("all_states",0)
		elif cmd.get("all_states")=='off':
			self.toggleStatesCaption='Toggle states: ON'
			cmd.set("all_states",1)
		self.cmd.refresh_wizard()
	
##########################################################################################
#superpose																				 #
##########################################################################################
	def superpose(self):
		#get the position of the selected residue's O atom
		stored.xyz = []
		if self.label.modifiedAA:
			cmd.iterate_state(1,self.residue1_name+" & name O","stored.xyz.append([x,y,z])")
		args=[]
		i = 0
		while i < len(self.label.atomsForSuperposition):
			args.append("%s & name %s" %(self.label.pymolName, self.label.atomsForSuperposition[i]))
			args.append("%s & name %s" %(self.residue1_name, self.label.atomsForSuperposition[i]))
			i+=1
		print(args)
		if cmd.pair_fit(*args):
			#set the label's O atom to the stored position
			if self.label.modifiedAA:
				cmd.alter_state(1,self.label.pymolName+" & name O","(x,y,z)=stored.xyz.pop(0)")
			return True
		else:
			return False

	def copyAndMove(self):
		newlabel=self.object_prefix+str(self.numberOfLabel)	   
		cmd.copy(newlabel+"_copied", self.picked_object1)
		self.superpose(newlabel+"_copied", self.residue2_name)
	
	def createSelectionMacro(self, selection):
		#selection+="& name CA"
		obj_list = cmd.get_names('objects')
		selectedObject=""
		for object in obj_list:
			if cmd.get_type(object)=="object:molecule":
				if cmd.count_atoms("(%s and (sele))"%(object)):
					#print self.picked_object1
					selectedObject=object
					break
		my_dict = { 'my_list' : [] }
		cmd.iterate(selection, "my_list.append((segi,chain,resi,resn))", space=my_dict)
		my_list = my_dict['my_list']
		#print my_list
		macro = "%s-%s-%s-%s-%s" % (selectedObject, my_list[0][0], my_list[0][1], my_list[0][2] ,my_list[0][3])
		return macro

	def createPseudoatom (self, coordinates, objectName, state):
		x=float(coordinates[0])
		y=float(coordinates[1])
		z=float(coordinates[2])
		posString="[%3.2f,%3.2f,%3.2f]" % (x,y,z)
		cmd.pseudoatom(pos=posString, object=objectName, state=state)	
	
##########################################################################################
#fastMtsslify																			 #
##########################################################################################
	def fastMtsslify(self, environmentAtoms):	#generate and check conformations semi-systematically
		#reference atoms are to detect internal clashes
		referenceAtoms=numpy.copy(self.label.movingAtoms)
		refDist=scipy.spatial.distance.cdist(referenceAtoms, referenceAtoms)
		
		#search settings for this label
		maxNtries=self.label.numberOfTries[self.thoroughness]
		numberToFind=self.label.numberToFind[self.thoroughness]
		found=0
		ntries=0
		axis = numpy.zeros(shape=(2,3)) 
		print("Trying to find conformations for label %s with vdW restraints: %s" % (self.label.pymolName, self.vdwRestraints))
		while found < numberToFind and ntries < maxNtries:
			self.label.movingAtoms = numpy.copy(referenceAtoms)
			if self.label.rotate:
				for chi in range (1, self.label.numberOfRotatingBonds+1):
					#print chi
					translationVector=numpy.array(self.label.movingAtoms[self.label.rotationInfo[str(chi)][0]])
					if chi == 1:
						backtranslationVector=translationVector
					else:
						backtranslationVector+=translationVector
					self.label.movingAtoms-=translationVector
					axis[0] = self.label.movingAtoms[self.label.rotationInfo[str(chi)][0]]
					axis[1] = self.label.movingAtoms[self.label.rotationInfo[str(chi)][0]+1]
					#rotate moving atoms around axis
					if not self.label.rotationInfo[str(chi)][2]:
						angle = generateRandomChiAngle()
					else:
						angle = generatePeptideChiAngle()
					rotationMatrix=setupRotationMatrix(angle, axis[1])
					self.label.movingAtoms[self.label.rotationInfo[str(chi)][1]]=rotatePoints(self.label.movingAtoms[self.label.rotationInfo[str(chi)][1]], rotationMatrix)
				self.label.movingAtoms+=backtranslationVector
			
			if not quickClash(self.label.movingAtoms[self.label.clashAtoms], environmentAtoms, self.cutoff, self.clashes) or not self.label.rotate:
				if not internalClash2(self.label.movingAtoms, refDist):
					found+=1
					print(found, end=' ')
					#only switch on snuggly fit search for "painstaking"
					if self.thoroughness == "painstaking":
						vdw = numberOfVdwContacts(self.label.movingAtoms[self.label.clashAtoms], environmentAtoms, self.cutoff)
						thisConformation = [found, vdw]
						self.conformationList.append(thisConformation)	
					self.createRotamer(found)
				else:
					print("i", end=' ')
			else:
				print(".", end=' ')
			ntries+=1
		results = [found, ntries]
		return results
	
##########################################################################################
#create Rotamer																				 #
##########################################################################################
	def createRotamer(self, found):
		for i in range (0, len(self.label.movingAtoms)):
			#print self.label.movingAtoms
			#print self.label.movingAtoms.tolist()
			stored.xyz = []
			if self.currentLabel == "Rx":
				stored.xyz = self.label.movingAtoms.tolist()
			else:
				stored.xyz = self.label.movingAtoms[i]
			#print self.label.pymolName
			#print stored.xyz
			try:
				cmd.alter_state(1,self.label.pymolName +"& name " + self.label.atomNames[i],"(x,y,z)=stored.xyz")
			except:
				pass
		cmd.create("%s_%s" %(self.residue1_name, self.label.identifier), self.label.pymolName, 1, found)

##########################################################################################
#create RxRotamer																				 #
##########################################################################################
	def createRxRotamer(self, solutions):
		objectName = "%s-%s_%s" %(self.residue1_name, self.residue2_name, self.label.identifier)
		for i in range (1, len(solutions)):
			self.createPseudoatom (solutions[i], objectName, i)


	
##########################################################################################
#create Snuggly fits																	 #
##########################################################################################
	def createSnugglyFitConformations(self): #select conformations with the best fit to the molecular surface, rank them and create an object
		print("Snuggliest fit(s):")
		#calculate average atom count of all conformations
		atomCountSum=0
		snugglyFitList=[]
		bestSnugglyFitAtomCount=0
		for x in range (0,len(self.conformationList)):
			thisConformation=self.conformationList[x]
			atomCountSum+=thisConformation[1]
			if thisConformation[1] > bestSnugglyFitAtomCount:
				bestSnugglyFitAtomCount=thisConformation[1]
		averageAtomCount=atomCountSum/len(self.conformationList)
		
		#generate snugglyFitList: take only those conformations whose atom count is > 0.75 of top peak and higher than the average
		counter=1
		for x in range (0,len(self.conformationList)):
			thisConformation=self.conformationList[x]
			if thisConformation[1] > 0.75 * bestSnugglyFitAtomCount and thisConformation[1] > averageAtomCount:
				snugglyFitList.append({'vdw':thisConformation[1],
									   'state':thisConformation[0]})
				#print "%s_snuggly_%s" %(self.residue1_name, self.label.identifier), "%s_%s" %(self.residue1_name, self.label.identifier), thisConformation[0], counter
				cmd.create("%s_snuggly_%s" %(self.residue1_name, self.label.identifier), "%s_%s" %(self.residue1_name, self.label.identifier), thisConformation[0], counter)
				counter+=1		
		#sort snugglyFitList so that best fitting conformation is on top   
		snugglyFitList = sorted(snugglyFitList, key=itemgetter('vdw'))
		snugglyFitList.reverse()
		#print out the list
		if len(snugglyFitList)>0:
			for i in range (0,len(snugglyFitList)):
				print("Conformation %i: \t\t%i \t\t\t vdW contacts" %(snugglyFitList[i]['state'],snugglyFitList[i]['vdw']))
		print("Average vdW contacts of all possible conformations: ",averageAtomCount)
	
	
	def calculateCone(self, ca, cb, environmentAtoms, numberOfAtoms=500):
		# generate umbrella of trial atoms
		#print self.label
		p = 2*self.label.trialAtomSphereRadius * numpy.random.rand(3, numberOfAtoms)- self.label.trialAtomSphereRadius
		p = p[:, sum(p* p, 0)** .5 <= self.label.trialAtomSphereRadius]
		p = p.transpose()
		p = p + cb
		distances = quick_map(ca, p)
		indices = numpy.where(numpy.any(distances < self.label.exclusionSphereRadius, axis=1))
		p = numpy.delete(p, indices, 0)
	
		#compute distances between trial sphere and environment
		distances = quick_map(p, environmentAtoms)
		#check for clashes and discard clashing atoms from sphere
		indices = numpy.where(numpy.any(distances < 3.5, axis=1))
		solutions = numpy.delete(p, indices,0)
		#print solutions
		numberOfConeAtoms = numpy.shape(solutions)[0]
	
		#plotting the spheres in Pymol
		ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
# 		for atom in p:
# 			self.createPseudoatom (atom, "cone"+ident, 1)		
# # 		
# 		for atom in solutions:
# 			self.createPseudoatom (atom, "solution"+ident, 1)
# 	
# 	
# 		averageConeCoordinate = numpy.average(solutions,axis=0)
# 		self.createPseudoatom(averageConeCoordinate, "avgCone"+ident, 1)
		if numberOfConeAtoms < 1:
			return ca, numberOfConeAtoms, solutions
		else:
			#calculate the average cone coordinate
			averageConeCoordinate = numpy.average(solutions,axis=0)
			return averageConeCoordinate, numberOfConeAtoms, solutions
	
##########################################################################################
#final cosmetics																		 #
##########################################################################################
	def finalCosmetics(self, found): #make everything look nice
		if found >= 1:
			print("Found %i conformations." %found)
			#show all conformations and do some coloring
			cmd.set("all_states",1)
			self.toggleStatesCaption='Toggle states: ON'
			cmd.color("blue","%s_%s" %(self.residue1_name, self.label.identifier))
			cmd.color("red", "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.highlight))
			util.cnc("%s_%s" %(self.residue1_name, self.label.identifier))
			cmd.disable(self.label.pymolName) 
			cmd.show("spheres", "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.highlight))
			cmd.set("sphere_scale", "0.2", "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.highlight))
			if self.currentLabel == "Rx":
				identifierLabel="%s-%s|%s|only N1 atoms!" %(self.residue1_name, self.residue2_name, self.currentLabel)
			else:
				identifierLabel="%s|%s|%s|%s" %(self.residue1_name, self.currentLabel, self.vdwRestraints, self.thoroughness)
			#print identifierLabel
			#mark average N1 position with pseudoatom
			stored.label = []
			cmd.iterate_state(0, "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.spinLocation), 'stored.label.append((x,y,z))')
			atoms1=numpy.array(stored.label)
			#create pseudoatom at average coordinate of each ensemble
			avgAtoms=numpy.average(atoms1,axis=0)
			self.createPseudoatom (avgAtoms, "%s_label" %(self.residue1_name), 1)
			cmd.set("sphere_scale", "0.5", "%s_label" %(self.residue1_name))
			cmd.label("%s_label" %(self.residue1_name), repr(identifierLabel))
			cmd.show("label")
			cmd.show("spheres", "name PS1")
			#print label+"*"+","+"labelEnvironment_"+label
			cmd.delete(self.label.pymolName+"*")
			cmd.group("%s%s" %(self.object_prefix, str(self.numberOfLabel)), "%s*, labelEnvironment_%s,%s*" %(self.label.pymolName, self.label.pymolName, self.residue1_name))
		else:
			print("Sorry, I did not find anything. Your options are:\n1) Try again, maybe with increased thoroughness,\n2) change vdW restraints to 'loose'")
			cmd.delete("%s*" %(self.label.pymolName))

##########################################################################################
#start wizard																			 #
##########################################################################################
def open_wizard():
	wiz = MtsslWizard()
	cmd.set_wizard(wiz)

##########################################################################################
#various functions																		 #
##########################################################################################

def quick_map(atoms1, atoms2):
	# if there is only one atom it has to be duplicated for quick_dist2 to work
	duplicated = False
	if len(numpy.shape(atoms1)) == 1:
		duplicated = True
		atoms1=numpy.tile(atoms1, (2,1))
	if len(numpy.shape(atoms2)) == 1:
		duplicated = True
		atoms2=numpy.tile(atoms2, (2,1))
	
	dist=scipy.spatial.distance.cdist(atoms1, atoms2)

	#remove the duplication depending on which selection contained the single atom
	if duplicated and dist.shape[0] == 2 and not dist.shape[1] == 2:
		dist=numpy.reshape(dist[0,:], (-1, 1))

	elif duplicated and dist.shape[1] == 2 and not dist.shape[0] == 2:
		dist=numpy.reshape(dist[:,0], (-1, 1))

	elif duplicated and dist.shape[0] == 2 and dist.shape[1] == 2:
		dist=numpy.reshape(dist[:1,0], (-1, 1))
	#else:
	#	dist=numpy.reshape(dist, (-1, 1))
	return dist #+6

def open_Plotter():
    # initialize window (roota)
    global plotWindow
    plotWindow = Tk()
    plotWindow.title(' Plotter by gha')
    global plotter
    plotter = Plotter(plotWindow)


def trialSphere():
	p= 2* numpy.random.rand(3, 1e4)- 1
	p= p[:, sum(p* p, 0)** .5<= 1]
	p=p*7
	sphere = numpy.hstack((p[0],p[1],p[2]))
	#print sphere
	return sphere
	

def generateRandomChiAngle():
	chi=random.random()*360.0
	return chi

def generatePeptideChiAngle():
	deltaChi=numpy.random.randint(-10,10)
	if random.choice([True, False]):
		chi=180+deltaChi
	else:
		chi=0+deltaChi
	return chi

def calculateStatistics2(distances):
	statisticsResult=""
	#statistics
	average = numpy.average(distances)
	median = numpy.median(distances)
	stddev = numpy.std(distances)
	longest = numpy.amax(distances)
	shortest = numpy.amin(distances)
	statisticsResult+= "Average of distribution: %3.2f\n" %average
	statisticsResult+= "Median of distribution: %3.2f\n" %median
	statisticsResult+= "Std. dev. of distribution: %3.2f\n" %stddev
	statisticsResult+= "Shortest distance: %3.2f\n" % shortest
	statisticsResult+= "Longest distance: %3.2f" %longest
	return statisticsResult

def numberOfVdwContacts(atoms, environmentAtoms, cutoff):
	dist=scipy.spatial.distance.cdist(environmentAtoms, atoms)
	vdwContacts=len(dist[numpy.nonzero((dist > cutoff) & (dist < 4.5 ))])
	return vdwContacts
	
def quick_dist2(atoms1, atoms2):
	# if there is only one atom it has to be duplicated for quick_dist2 to work
	duplicated = False
	if len(atoms1) == 1:
		#print "A"
		duplicated = True
		atoms1=numpy.tile(atoms1, (2,1))
	if len(atoms2) == 1:
		#print "B"
		duplicated = True
		atoms2=numpy.tile(atoms2, (2,1))
	#print atoms1
	#print atoms2	
	#take random sample if too many atoms
	if len(atoms1) > 250:
		#atoms1=atoms1[numpy.random.randint(atoms1.shape[0], size = 10),:]
		atoms1=atoms1[numpy.random.permutation(atoms1.shape[0])[:250]]
	if len(atoms2) > 250:
		#atoms2=atoms2[numpy.random.randint(atoms2.shape[0], size = 10),:]
		atoms2=atoms2[numpy.random.permutation(atoms2.shape[0])[:250]]
	dist=scipy.spatial.distance.cdist(atoms1, atoms2)
	
	#remove the duplication depending on which selection contained the single atom
	if duplicated and dist.shape[0] == 2 and not dist.shape[1] == 2:
		dist=numpy.reshape(dist[0,:], (-1, 1))
	
	elif duplicated and dist.shape[1] == 2 and not dist.shape[0] == 2:
		dist=numpy.reshape(dist[:,0], (-1, 1))
	
	elif duplicated and dist.shape[0] == 2 and dist.shape[1] == 2:
		dist=numpy.reshape(dist[:1,0], (-1, 1))
	else:
		dist=numpy.reshape(dist, (-1, 1))
	
	return dist

def internalClash2(atoms, refDist):
	#distances in new rotamer
	dist=scipy.spatial.distance.cdist(atoms, atoms)
	#create Boolean array with elements that describe if a distance changes or not
	changingDistances = numpy.absolute(numpy.round(numpy.subtract(dist,refDist),2)) > 0
	#multiply by Boolean array to make all constant distances zero
	dist=changingDistances*dist
	#check for internal clashes
	internalClashes=dist[numpy.nonzero((dist < internalClash_cutoff) & (dist > 0))]
	#print internalClashes
	if len(internalClashes) > 0:
		return True
	else:
		return False

def quickClash(rotatedAtoms, environmentAtoms, cutoff, maxClash):
	dist=scipy.spatial.distance.cdist(environmentAtoms, rotatedAtoms)
	clashes=len(numpy.nonzero(dist < cutoff)[0])
	if clashes > maxClash:
		return True
	else:
		return False
	return False

def rotatePoints(points, rotationMatrix):
	rotatedPoints=[]
	#print points
	#print rotationMatrix
	for point in points:
		#add 1 for multiplication with 4x4 matrix
		point = numpy.append(point, 1)
		rotatedPoint=numpy.dot(rotationMatrix,point)
		#remove 1 again
		rotatedPoint=numpy.delete(rotatedPoint, 3)
		rotatedPoints.append(rotatedPoint)
	return numpy.array(rotatedPoints)
		
def setupRotationMatrix(angle, axisPoint):
	#print axisPoint
	u = axisPoint[0]
	v = axisPoint[1]
	w = axisPoint[2]
	L = (u*u + v * v + w * w)
	angle = angle * numpy.pi / 180.0
	u2 = u * u
	v2 = v * v
	w2 = w * w
	rotationMatrix = numpy.zeros((4,4))
	rotationMatrix[0][0] = (u2 + (v2 + w2) * numpy.cos(angle)) / L
	rotationMatrix[0][1] = (u * v * (1 - numpy.cos(angle)) - w * numpy.sqrt(L) * numpy.sin(angle)) / L
	rotationMatrix[0][2] = (u * w * (1 - numpy.cos(angle)) + v * numpy.sqrt(L) * numpy.sin(angle)) / L
	rotationMatrix[0][3] = 0.0
	
	rotationMatrix[1][0] = (u * v * (1 - numpy.cos(angle)) + w * numpy.sqrt(L) * numpy.sin(angle)) / L
	rotationMatrix[1][1] = (v2 + (u2 + w2) * numpy.cos(angle)) / L
	rotationMatrix[1][2] = (v * w * (1 - numpy.cos(angle)) - u * numpy.sqrt(L) * numpy.sin(angle)) / L
	rotationMatrix[1][3] = 0.0
	
	rotationMatrix[2][0] = (u * w * (1 - numpy.cos(angle)) - v * numpy.sqrt(L) * numpy.sin(angle)) / L
	rotationMatrix[2][1] = (v * w * (1 - numpy.cos(angle)) + u * numpy.sqrt(L) * numpy.sin(angle)) / L
	rotationMatrix[2][2] = (w2 + (u2 + v2) * numpy.cos(angle)) / L
	rotationMatrix[2][3] = 0.0
	
	rotationMatrix[3][0] = 0.0
	rotationMatrix[3][1] = 0.0
	rotationMatrix[3][2] = 0.0
	rotationMatrix[3][3] = 1.0
	return rotationMatrix

#########################################################################################
#Label classes           																#
#after adding a new class, don't forget to make changes in 'run',and add it to the GUI! # 
#########################################################################################
class MtsslLabel:
	identifier = "M-T-S-S-L"
	modifiedAA = True
	numberOfRotatingBonds = 5
	numberOfAtoms = 18
	rotate = True
	rotationInfo = {'1': [3,slice(5, numberOfAtoms), False], '2': [4,slice(6, numberOfAtoms), False], '3': [5,slice(7, numberOfAtoms), False], '4': [6,slice(8, numberOfAtoms), False], '5': [7,slice(9, numberOfAtoms), False]}
	clashAtoms = slice(5, numberOfAtoms)
	radius = 13.0
	atomNames =	 ['N', 'O', 'C', 'CA', 'CB', 'SG', 'SD', 'CE', 'C3', 'O1', 'C2', 'N1', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
				#  0	1	 2	  3		4	  5		6	  7		8	  9		10	  11	12	  13	14	  15	16	  17 
	
	unsortedAtomNames = ['N', 'CA', 'C', 'O', 'CB', 'SG', 'SD', 'CE', 'N1', 'O1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
	movingAtoms=[]
	spinLocation = 'N1'
	highlight = 'O1'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 2000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = ""
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    MTSSL\n                
ATOM      1 N    R1A A   1       0.201  -0.038  -0.149  1.00 20.00           N\n
ATOM      2 CA   R1A A   1       1.258   1.007  -0.271  1.00 20.00           C\n
ATOM      4 CB   R1A A   1       2.056   0.796  -1.554  1.00 20.00           C\n
ATOM      7 SG   R1A A   1       3.667   1.492  -1.387  1.00 20.00           S\n
ATOM      8 SD   R1A A   1       4.546   1.587  -3.180  1.00 20.00           S\n
ATOM      9 CE   R1A A   1       5.573   3.020  -3.244  1.00 20.00           C\n
ATOM     12 C3   R1A A   1       6.644   3.007  -4.321  1.00 20.00           C\n
ATOM     13 C4   R1A A   1       7.349   4.109  -4.593  1.00 20.00           C\n
ATOM     15 C5   R1A A   1       8.361   3.885  -5.680  1.00 20.00           C\n
ATOM     16 C7   R1A A   1       9.758   4.118  -5.119  1.00 20.00           C\n
ATOM     20 C6   R1A A   1       8.092   4.808  -6.859  1.00 20.00           C\n
ATOM     24 N1   R1A A   1       8.144   2.482  -5.989  1.00 20.00           N\n
ATOM     25 O1   R1A A   1       8.792   1.889  -6.835  1.00 20.00           O\n
ATOM     26 C2   R1A A   1       7.092   1.857  -5.197  1.00 20.00           C\n
ATOM     27 C8   R1A A   1       7.642   0.712  -4.360  1.00 20.00           C\n
ATOM     31 C9   R1A A   1       5.961   1.403  -6.123  1.00 20.00           C\n
ATOM     35 C    R1A A   1       0.670   2.388  -0.261  1.00 20.00           C\n
ATOM     36 O    R1A A   1      -0.298   2.670  -0.967  1.00 20.00           O\n"""
	
	def __init__(self, pymolName):
		self.pymolName = pymolName

	
class ProxylLabel:
	identifier = "P-R-O-X-Y-L"
	modifiedAA = True
	numberOfRotatingBonds = 6
	numberOfAtoms = 20
	rotate = True
	rotationInfo = {'1': [3,slice(5, numberOfAtoms), False], '2': [4,slice(6, numberOfAtoms), False], '3': [5,slice(7, numberOfAtoms), False], '4': [6,slice(8, numberOfAtoms), False], '5': [7,slice(9, numberOfAtoms-1), False], '6': [8,slice(9, numberOfAtoms-1), False]}
	clashAtoms = slice(5, numberOfAtoms)
	radius = 13.0
	atomNames =	 ['N', 'O', 'C', 'CA', 'CB', 'SG', 'C1', 'C2', 'N2', 'C3', 'N1', 'O1', 'C10','C4', 'C5', 'C6', 'C7', 'C8', 'C9',  'O3']
				#  0	1	 2	  3		4	  5		6	  7		8	  9		10	  11	12	  13	14	  15	16	  17	18	   19
	movingAtoms=[]
	spinLocation = 'N1'
	highlight = 'O1'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 20000, 'quick search': 3000}
	info = ""
	errorMessage = ""
	trialAtomSphereRadius = 9.5
	exclusionSphereRadius = 9.5
	pdbStr = """HEADER    PROXYL\n
ATOM      1 N    IA1 A   1       0.727  -0.687   0.805  1.00 20.00           N\n
ATOM      2 CA   IA1 A   1      -0.684  -0.504   1.000  1.00 20.00           C\n
ATOM      4 C    IA1 A   1      -1.116   0.797   0.371  1.00 20.00           C\n
ATOM      7 O    IA1 A   1      -2.273   1.140   0.424  1.00 20.00           O\n
ATOM      8 CB   IA1 A   1      -0.990  -0.471   2.492  1.00 20.00           C\n
ATOM      9 SG   IA1 A   1      -2.133  -1.819   2.906  1.00 20.00           S\n
ATOM     12 C1   IA1 A   1      -2.688  -2.044   4.619  1.00 20.00           C\n
ATOM     13 C2   IA1 A   1      -3.625  -3.224   4.694  1.00 20.00           C\n
ATOM     15 O3   IA1 A   1      -3.891  -3.849   3.694  1.00 20.00           O\n
ATOM     16 N2   IA1 A   1      -4.175  -3.587   5.894  1.00 20.00           N\n
ATOM     20 C3   IA1 A   1      -5.075  -4.721   5.966  1.00 20.00           C\n
ATOM     24 C4   IA1 A   1      -4.559  -5.760   6.953  1.00 20.00           C\n
ATOM     25 C5   IA1 A   1      -5.811  -6.228   7.683  1.00 20.00           C\n
ATOM     26 N1   IA1 A   1      -6.929  -5.507   7.146  1.00 20.00           N\n
ATOM     27 C6   IA1 A   1      -6.424  -4.306   6.540  1.00 20.00           C\n
ATOM     31 C7   IA1 A   1      -6.006  -7.727   7.496  1.00 20.00           C\n
ATOM     35 C8   IA1 A   1      -5.678  -5.942   9.173  1.00 20.00           C\n
ATOM     36 O1   IA1 A   1      -8.190  -5.879   7.200  1.00 20.00           O\n
ATOM     37 C9   IA1 A   1      -6.309  -3.153   7.529  1.00 20.00           C\n
ATOM     38 C10  IA1 A   1      -7.356  -3.824   5.436  1.00 20.00           C\n"""
	
	def __init__(self, pymolName):
		self.pymolName = pymolName

		
class UripLabel:
	identifier = "U-R-I-P"
	modifiedAA = False
	numberOfRotatingBonds = 4
	numberOfAtoms = 18
	rotate = True
	rotationInfo = {'1': [2,slice(4, numberOfAtoms), False], '2': [3,slice(5, numberOfAtoms), False],'3': [4,slice(6, numberOfAtoms-1), False],'4': [5,slice(7, numberOfAtoms-1), False]}
	clashAtoms = slice(4, numberOfAtoms)
	radius = 13.0
	atomNames =	 ["O4'", "C3'", "C2'", 'NS1', 'CS1', 'NS2', 'CS7', 'CS2', 'CS10', 'CS11',  'CS3', 'CS5', 'CS6',  'CS8', 'CS9',   'NS3', 'OS1', 'OS2']
	movingAtoms=[]
	spinLocation = 'NS3'
	highlight = 'OS1'
	atomsForSuperposition = ["C2'","C3'","O4'"]
	defaultVdw = "loose"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 20000, 'quick search': 3000}
	info = "\nFor the Urip label, the vdW restraints are by default set to 'loose' to account for possible polar interactions between the amide bonds and the DNA backbone."
	errorMessage = "Check atom nomenclature. The ribose atoms are sometimes called C2* instead of C2'\nThis can be changed e.g. by 'alter all, name=name.replace(\"*\",\"'\")' in PyMOL."
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    URIPSL\n
ATOM      2  O4' URI A   6     -22.786  68.384   8.719  1.00  0.00           O  
ATOM      3  C3' URI A   6     -23.259  70.361   9.877  1.00  0.00           C  
ATOM      5  C2' URI A   6     -22.457  69.433  10.789  1.00  0.00           C  
ATOM      7  CS1 URI A   6     -23.026  70.760  12.707  0.00  0.00           C  
ATOM      8  CS7 URI A   6     -23.409  72.326  14.678  0.00  0.00           C  
ATOM      9  NS1 URI A   6     -23.018  69.401  12.147  1.00  0.00           N  
ATOM     10  OS1 URI A   6     -22.138  73.362  17.681  0.00  0.00           O  
ATOM     11  CS2 URI A   6     -24.196  72.339  15.999  0.00  0.00           C  
ATOM     12  NS2 URI A   6     -23.400  70.968  14.114  0.00  0.00           N  
ATOM     13  OS2 URI A   6     -22.728  71.707  12.019  0.00  0.00           O  
ATOM     14  CS3 URI A   6     -24.179  73.739  16.634  0.00  0.00           C  
ATOM     15  NS3 URI A   6     -22.776  74.200  16.758  0.00  0.00           N  
ATOM     16  CS5 URI A   6     -21.952  74.215  15.527  0.00  0.00           C  
ATOM     17  CS6 URI A   6     -21.968  72.816  14.891  0.00  0.00           C
ATOM      4  CS8 URI A   6     -22.528  75.240  14.533  0.00  0.00           C  
ATOM      5  CS9 URI A   6     -20.505  74.607  15.878  0.00  0.00           C  
ATOM      6 CS10 URI A   6     -24.969  74.719  15.747  0.00  0.00           C  
ATOM      7 CS11 URI A   6     -24.829  73.681  18.029  0.00  0.00           C\n"""
	
	def __init__(self, pymolName):
		self.pymolName = pymolName
	
class CLabel:
	identifier = "C-L-A-B-E-L"
	modifiedAA = False
	numberOfRotatingBonds = 0
	numberOfAtoms = 23
	rotate = False
	rotationInfo = {}
	clashAtoms = slice(4, numberOfAtoms)
	radius = 20
	atomNames =	 ['C1', 'N1', 'C2', 'O2', 'N3', 'C5', 'C6', 'N7', 'C8', 'C9', 'O10', 'C11', 'C12', 'C13', 'C14', 'C15', 'N16', 'C17', 'C19', 'C20', 'O21', 'C22', 'C23']
	movingAtoms=[]
	spinLocation = 'N16'
	highlight = 'O21'
	atomsForSuperposition = ["N1", "C2", "O2", "N3"]
	#atomsForSuperposition = ["C2'","C3'","O4'"]
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1, 'thorough search': 1, 'quick search': 1}
	numberOfTries = {'painstaking': 1, 'thorough search': 1, 'quick search': 1}
	info = "\nThe C label is only superimposed onto e.g. dC but does not move.\n"
	errorMessage = "This label does not superpose onto A or G!"
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    CLABEL\n
HETATM   10  C1  EXC B   2       8.773   1.726  20.049  1.00 15.03           C  
HETATM   11  N3  EXC B   2       8.201   2.413  21.056  1.00 12.25           N  
HETATM   12  C2  EXC B   2       7.676   1.752  22.102  1.00 10.61           C  
HETATM   13  N1  EXC B   2       7.709   0.407  22.179  1.00 11.48           N  
HETATM   14  C5  EXC B   2       8.265  -0.311  21.201  1.00 11.57           C  
HETATM   15  C6  EXC B   2       8.806   0.340  20.112  1.00 16.49           C  
HETATM   16  N7  EXC B   2       9.308   2.376  18.989  1.00 12.49           N  
HETATM   17  C8  EXC B   2       9.883   1.677  17.984  1.00 16.11           C  
HETATM   18  C9  EXC B   2       9.920   0.280  18.035  1.00 16.30           C  
HETATM   19  O10 EXC B   2       9.375  -0.375  19.097  1.00 15.90           O  
HETATM   20  C11 EXC B   2      10.514  -0.465  17.020  1.00 18.81           C  
HETATM   21  C12 EXC B   2      11.047   0.229  15.945  1.00 20.67           C  
HETATM   22  C13 EXC B   2      11.010   1.597  15.888  1.00 24.21           C  
HETATM   23  C14 EXC B   2      10.433   2.357  16.897  1.00 19.54           C  
HETATM   24  C15 EXC B   2      11.663   2.110  14.644  1.00 27.86           C  
HETATM   25  N16 EXC B   2      12.121   0.858  14.038  1.00 29.97           N  
HETATM   26  C17 EXC B   2      11.734  -0.359  14.747  1.00 23.17           C  
HETATM   27  O2  EXC B   2       7.152   2.387  23.031  1.00 11.16           O  
HETATM   28  C19 EXC B   2      10.648   2.835  13.770  1.00 21.72           C  
HETATM   29  C20 EXC B   2      12.837   3.004  15.031  1.00 21.39           C  
HETATM   30  O21 EXC B   2      12.786   0.833  13.018  1.00 39.82           O  
HETATM   31  C22 EXC B   2      10.769  -1.222  13.943  1.00 31.01           C  
HETATM   32  C23 EXC B   2      12.956  -1.146  15.215  1.00 19.65           C  \n"""
	
	def __init__(self, pymolName):
		self.pymolName = pymolName

class Dota1Label:
	identifier = "D-O-T-A-1"
	modifiedAA = True
	numberOfRotatingBonds = 6
	numberOfAtoms = 38
	rotate = True
	rotationInfo = {'1': [3,slice(5, numberOfAtoms), False], '2': [4,slice(6, numberOfAtoms), False], '3': [5,slice(7, numberOfAtoms), False], '4': [6,slice(8, numberOfAtoms), False], '5': [7,slice(9, numberOfAtoms), False], '6': [8,slice(10, numberOfAtoms), False]}
	clashAtoms = slice(5, numberOfAtoms)
	radius = 18.0
	atomNames =	 ['N', 'O',  'C', 'CA', 'CB',  'SG', 'SD', 'C17','C16','N5', 'C15', 'C3', 'N3', 'O3', 'C4', 'N4', 'O4', 'C5', 'N2', 'O5', 'C6', 'O6', 'C7', 'O7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'O2',  'C2',  'O1',  'C18', 'Gd', 'C1','N1']
				#  0	1	  2	   3	 4	    5	  6	    7	  8	    9	  10	 11    12	 13    14	 15    16	 17    18    19    20    21    22    23    24    25    26     27     28     29     30     31     32     33     34     35    36   37 
	
	movingAtoms=[]
	spinLocation = 'Gd'
	highlight = 'Gd'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = ""
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    DTA\n
ATOM      1 N    DTA A   1       0.201  -0.038  -0.149  1.00 20.00           N 
ATOM      2 C    R1A A   1       0.670   2.388  -0.261  1.00 20.00           C
ATOM      3 O    R1A A   1      -0.298   2.670  -0.967  1.00 20.00           O
ATOM      4 CA   DTA A   1       1.258   1.007  -0.271  1.00 20.00           C 
ATOM      5 CB   DTA A   1       2.056   0.796  -1.554  1.00 20.00           C 
ATOM      6 SG   DTA A   1       3.667   1.492  -1.387  1.00 20.00           S 
ATOM      7 SD   DTA A   1       4.546   1.587  -3.180  1.00 20.00           S 
ATOM      8 N1   DTA A   1      10.664   6.956   2.914  0.00  0.00           N  
ATOM      9 N2   DTA A   1       8.088   7.390   1.217  0.00  0.00           N  
ATOM     10 N3   DTA A   1       8.072   4.476   0.262  0.00  0.00           N  
ATOM     11 N4   DTA A   1      10.642   3.874   2.091  0.00  0.00           N  
ATOM     12 N5   DTA A   1       7.540   4.129  -2.243  0.00  0.00           N  
ATOM     13 C1   DTA A   1       8.581   8.583   2.347  0.00  0.00           C  
ATOM     14 C2   DTA A   1       9.535   8.072   3.473  0.00  0.00           C  
ATOM     15 C3   DTA A   1       6.635   6.520   1.215  0.00  0.00           C  
ATOM     16 C4   DTA A   1       6.735   5.019   0.751  0.00  0.00           C  
ATOM     17 C5   DTA A   1      11.218   5.797   4.032  0.00  0.00           C  
ATOM     18 C6   DTA A   1      11.595   4.402   3.412  0.00  0.00           C  
ATOM     19 C7   DTA A   1       9.321   2.786   1.945  0.00  0.00           C  
ATOM     20 C8   DTA A   1       8.376   3.044   0.699  0.00  0.00           C  
ATOM     21 C9   DTA A   1      11.962   7.867   2.554  0.00  0.00           C  
ATOM     22 C10  DTA A   1      12.595   7.180   1.363  0.00  0.00           C  
ATOM     23 C11  DTA A   1      11.770   3.656  -0.089  0.00  0.00           C  
ATOM     24 C12  DTA A   1      11.741   3.037   1.290  0.00  0.00           C  
ATOM     25 C13  DTA A   1       7.919   8.327  -0.072  0.00  0.00           C  
ATOM     26 C14  DTA A   1       9.320   8.374  -0.663  0.00  0.00           C  
ATOM     27 C15  DTA A   1       8.129   4.708  -1.156  0.00  0.00           C  
ATOM     28 C16  DTA A   1       6.694   2.935  -2.181  0.00  0.00           C  
ATOM     29 C17  DTA A   1       5.585   3.049  -3.239  0.00  0.00           C  
ATOM     30 C18  DTA A   1       8.591   5.328   2.867  0.00  0.00           C  
ATOM     31 O1   DTA A   1      12.279   6.121   0.827  0.00  0.00           O  
ATOM     32 O2   DTA A   1      13.723   7.836   0.992  0.00  0.00           O  
ATOM     33 O3   DTA A   1      10.918   4.338  -0.654  0.00  0.00           O  
ATOM     34 O4   DTA A   1      12.860   3.224  -0.772  0.00  0.00           O  
ATOM     35 O5   DTA A   1      10.366   7.900  -0.217  0.00  0.00           O  
ATOM     36 O6   DTA A   1       9.331   9.175  -1.757  0.00  0.00           O  
ATOM     37 O7   DTA A   1       8.864   5.668  -1.199  0.00  0.00           O  
ATOM     38 Gd   DTA A   1       9.943   5.855   0.862  0.00  0.00          GD\n"""
	
	def __init__(self, pymolName):
		self.pymolName = pymolName
	
class ByspLabel:
	identifier = "B-Y-S-P"
	modifiedAA = True
	numberOfRotatingBonds = 9
	numberOfAtoms = 23
	rotate = True
	rotationInfo = {'1': [3,slice(5, numberOfAtoms), False], '2': [4,slice(6, numberOfAtoms), False], '3': [5,slice(7, numberOfAtoms), False], '4': [6,slice(8, numberOfAtoms), False], '5': [7,slice(9, numberOfAtoms), False], '6': [17,slice(19, numberOfAtoms), False], '7': [18,slice(20, numberOfAtoms), False], '8': [19,slice(21, numberOfAtoms), False], '9': [20,slice(22, numberOfAtoms), False]}
	clashAtoms = slice(5, 5)
	radius = 16.0
	atomNames =	['N', 'O', 'C', 'CA', 'CB', 'SG1', 'SG2', 'CL1', 'CS3', 'CS2', 'CS8', 'CS9', 'NS1', 'OS1', 'CS5', 'CS6', 'CS7', 'CS4', 'CL2', 'SG3', 'SG4', 'CB2', 'CA1'] 
	#             0    1    2    3     4     5      6      7      8      9      10     11     12     13     14     15     16     17     18     19     20     21     22 
	
	movingAtoms=[]
	spinLocation = 'NS1'
	highlight = 'OS1'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = ""
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    BYSP\n                
ATOM     50  N   BYSP   16       2.664  32.900  33.516  1.00 37.45      A\n
ATOM     51  CA  BYSP   16       3.870  32.287  34.069  1.00 37.25      A\n
ATOM     52  CB  BYSP   16       3.521  31.459  35.330  1.00 41.00      A\n
ATOM     53  SG1 BYSP   16       2.337  30.121  34.923  1.00  0.00      A\n
ATOM     54  SG2 BYSP   16       3.362  28.406  35.157  1.00  0.00      A\n
ATOM     55  SG3 BYSP   16       7.434  28.620  30.786  1.00  0.00      A\n
ATOM     56  CL1 BYSP   16       4.846  28.679  34.193  1.00  0.00      A\n
ATOM     57  CS3 BYSP   16       5.725  27.479  33.984  1.00  0.00      A\n
ATOM     58  CS2 BYSP   16       6.175  26.549  35.100  1.00  0.00      A\n
ATOM     59  CS8 BYSP   16       6.974  27.286  36.213  1.00  0.00      A\n
ATOM     60  CS9 BYSP   16       4.967  25.810  35.626  1.00  0.00      A\n
ATOM     61  NS1 BYSP   16       7.009  25.616  34.335  1.00  0.00      A\n
ATOM     62  OS1 BYSP   16       7.555  24.642  34.879  1.00  0.00      A\n
ATOM     63  CS5 BYSP   16       7.106  25.902  32.887  1.00  0.00      A\n
ATOM     64  CS7 BYSP   16       8.512  26.343  32.589  1.00  0.00      A\n
ATOM     65  CS6 BYSP   16       6.715  24.752  32.010  1.00  0.00      A\n
ATOM     66  CS4 BYSP   16       6.225  27.124  32.813  1.00  0.00      A\n
ATOM     67  CL2 BYSP   16       5.969  27.823  31.527  1.00  0.00      A\n
ATOM     68  C   BYSP   16       4.541  31.404  33.041  1.00 35.65      A\n
ATOM     69  O   BYSP   16       5.642  30.960  33.232  1.00 35.84      A\n
ATOM     88  CA1 BYSP   16       9.638  32.113  31.010  1.00 32.61      A\n
ATOM     89  CB2 BYSP   16       9.444  30.675  31.495  1.00 35.74      A\n
ATOM     90  SG4 BYSP   16       7.715  30.249  31.870  1.00  0.00      A\n"""
	
	def __init__(self, pymolName):
		self.pymolName = pymolName
		
class PAcPheLabel:
	identifier = "p-A-c-P-h-e"
	modifiedAA = True
	numberOfRotatingBonds = 6
	numberOfAtoms = 26
	rotate = True
	rotationInfo = {'1': [3,slice(5, numberOfAtoms), False], '2': [4,slice(6, numberOfAtoms), False], '3': [10,slice(12, numberOfAtoms), True], '4': [12,slice(14, numberOfAtoms-1), False], '5': [13,slice(15, numberOfAtoms-1), False], '6': [14,slice(16, numberOfAtoms-1), False]}
	clashAtoms = slice(5, numberOfAtoms)
	radius = 15.0
	atomNames =	['N', 'O', 'C', 'CA', 'CB', 'C19', 'C18', 'C17', 'C25', 'C26', 'C16', 'C02', 'N03', 'O04', 'CE', 'C3', 'O1', 'C2', 'N1', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C01']
	           #  0	   1	2	  3		4	   5	  6      7      8      9     10     11     12     13     14     15    16    17    18    19    20    21    22    23    24    25
	
	movingAtoms=[]
	spinLocation = 'N1'
	highlight = 'O1'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = ""
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    LIG\n                
ATOM      1  N   LIG A   1       9.798  -1.128   4.576  1.00 20.00           N\n  
ATOM      2  CA  LIG A   1      10.446  -0.181   3.713  1.00 20.00           C\n 
ATOM      3  C   LIG A   1      11.904  -0.389   3.763  1.00 20.00           C\n  
ATOM      4  O   LIG A   1      12.516  -0.159   4.800  1.00 20.00           O\n  
ATOM      5  CB  LIG A   1      10.060   1.237   4.126  1.00 20.00           C\n  
ATOM      6  CE  LIG A   1       3.592   4.420  -0.998  1.00 20.00           C\n  
ATOM      7  C01 LIG A   1       7.779   4.832  -1.222  1.00 20.00           C\n  
ATOM      8  N1  LIG A   1       2.136   1.331  -2.662  1.00 20.00           N\n  
ATOM      9  O1  LIG A   1       1.410   0.363  -2.819  1.00 20.00           O\n  
ATOM     10  C02 LIG A   1       7.034   4.049  -0.163  1.00 20.00           C\n  
ATOM     11  C2  LIG A   1       2.033   2.237  -1.525  1.00 20.00           C\n  
ATOM     12  C3  LIG A   1       3.147   3.219  -1.815  1.00 20.00           C\n  
ATOM     13  N03 LIG A   1       5.702   3.974  -0.203  1.00 20.00           N\n  
ATOM     14  C4  LIG A   1       3.774   2.900  -2.952  1.00 20.00           C\n  
ATOM     15  O04 LIG A   1       5.013   4.622  -1.186  1.00 20.00           O\n  
ATOM     16  C5  LIG A   1       3.189   1.678  -3.600  1.00 20.00           C\n  
ATOM     17  C6  LIG A   1       2.669   2.009  -4.991  1.00 20.00           C\n  
ATOM     18  C7  LIG A   1       4.252   0.589  -3.666  1.00 20.00           C\n  
ATOM     19  C8  LIG A   1       2.276   1.503  -0.215  1.00 20.00           C\n  
ATOM     20  C9  LIG A   1       0.670   2.932  -1.552  1.00 20.00           C\n  
ATOM     21  C16 LIG A   1       7.811   3.339   0.966  1.00 20.00           C\n  
ATOM     22  C17 LIG A   1       7.130   2.612   1.931  1.00 20.00           C\n  
ATOM     23  C18 LIG A   1       7.853   1.931   2.956  1.00 20.00           C\n  
ATOM     24  C19 LIG A   1       9.254   1.997   2.987  1.00 20.00           C\n  
ATOM     25  C25 LIG A   1       9.939   2.726   2.021  1.00 20.00           C\n  
ATOM     26  C26 LIG A   1       9.217   3.406   0.998  1.00 20.00           C\n"""

	def __init__(self, pymolName):
		self.pymolName = pymolName



class TPA:
	identifier = "T-P-A"
	modifiedAA = False
	numberOfRotatingBonds = 1
	numberOfAtoms = 16
	rotate = True
	rotationInfo = {'1': [5,slice(7, numberOfAtoms), False]}
	clashAtoms = slice(8, numberOfAtoms)
	radius = 15.0
	atomNames =	['C4', 'C5', 'C6', 'O14', 'C13', 'C12', 'C15', 'C14', 'C7', 'C8', 'C9', 'O15', 'N1', 'C16', 'C10', 'C11']
	           #  0     1     2     3      4      5      6      7      8      9    10    11     12    13     14     15
	movingAtoms=[]
	spinLocation = 'N1'
	highlight = 'O15'
	atomsForSuperposition = ['C4','C5','C6']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = "Only works for Uridine, Cytosine or Thymine"
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    LIG\n                
ATOM      1  N1          1      12.207   2.205  -0.447  0.00  0.00           N  \n
ATOM      2  C4          1       6.030  -0.422   0.193  0.00  0.00           C  \n
ATOM      3  C5          1       6.969  -0.531   1.139  0.00  0.00           C  \n
ATOM      4  C6          1       6.716  -1.293   2.412  0.00  0.00           C  \n
ATOM      5  C14         1      11.211   1.509   1.573  0.00  0.00           C  \n
ATOM      6  C15         1      10.384   1.136   0.597  0.00  0.00           C  \n
ATOM      7  C16         1      10.968   1.457  -0.729  0.00  0.00           C  \n
ATOM      8  C7          1      12.440   2.119   1.006  0.00  0.00           C  \n
ATOM      9  C8          1      12.661   3.508   1.592  0.00  0.00           C  \n
ATOM     10  C9          1      13.679   1.295   1.330  0.00  0.00           C  \n
ATOM     11  C10         1      10.004   2.312  -1.541  0.00  0.00           C  \n
ATOM     12  C11         1      11.242   0.199  -1.542  0.00  0.00           C  \n
ATOM     13  C12         1       9.216   0.566   0.782  0.00  0.00           C  \n
ATOM     14  C13         1       8.137   0.040   0.954  0.00  0.00           C  \n
ATOM     15  O14         1       7.407  -1.167   3.394  0.00  0.00           O  \n
ATOM     16  O15         1      12.961   2.831  -1.326  0.00  0.00           O  \n"""

	def __init__(self, pymolName):
		self.pymolName = pymolName

class KznLabel:
	identifier = "K-Z-N"
	modifiedAA = False
	numberOfRotatingBonds = 2
	numberOfAtoms = 26
	rotate = True
	rotationInfo = {'1': [6,slice(8, numberOfAtoms), False], '2': [12,slice(14, numberOfAtoms), False]}
	clashAtoms = slice(8, numberOfAtoms)
	radius = 15.0
	atomNames =['N32','O34', 'O35', 'C4', 'C6', 'N01', 'C5', 'C12', 'C13', 'N29', 'N30', 'N14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'N27', 'O28', 'C33', ]
	           # 0     1     2     3      4      5      6      7      8      9      10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26	
	movingAtoms=[]
	spinLocation = 'N27'
	highlight = 'O28'
	atomsForSuperposition = ['C6','C5','C4']
	defaultVdw = "loose"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = "Only works for Uridine, Cytosine or Thymine"
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    KZN\n                
HETATM    1  C4 K ZN     1      -2.543  -2.393  -1.845  1.00 20.00           C  \n
HETATM    2  C5 K ZN     1      -2.032  -1.830  -0.590  1.00 20.00           C  \n
HETATM    3  C6 K ZN     1      -2.856  -1.892   0.628  1.00 20.00           C  \n
HETATM    4  N01 KZN A   1      -4.057  -2.654   0.634  1.00 20.00      A    N  \n
HETATM    5  C12 KZN A   1      -0.608  -1.301  -0.518  1.00 20.00      A    C  \n
HETATM    6  C13 KZN A   1      -0.079  -0.141  -1.166  1.00 20.00      A    C  \n
HETATM    7  N14 KZN A   1       1.224  -0.091  -0.891  1.00 20.00      A    N  \n
HETATM    8  C15 KZN A   1       2.123   0.944  -1.296  1.00 20.00      A    C  \n
HETATM    9  C16 KZN A   1       2.704   0.888  -2.525  1.00 20.00      A    C  \n
HETATM   10  C17 KZN A   1       3.734   1.751  -2.841  1.00 20.00      A    C  \n
HETATM   11  C18 KZN A   1       4.192   2.709  -1.873  1.00 20.00      A    C  \n
HETATM   12  C19 KZN A   1       5.607   3.224  -1.746  1.00 20.00      A    C  \n
HETATM   13  C20 KZN A   1       6.609   2.101  -2.054  1.00 20.00      A    C  \n
HETATM   14  C21 KZN A   1       5.827   4.399  -2.694  1.00 20.00      A    C  \n
HETATM   15  C22 KZN A   1       2.567   1.867  -0.368  1.00 20.00      A    C  \n
HETATM   16  C23 KZN A   1       3.639   2.764  -0.701  1.00 20.00      A    C  \n
HETATM   17  C24 KZN A   1       4.455   3.671   0.148  1.00 20.00      A    C  \n
HETATM   18  C25 KZN A   1       3.886   5.093   0.101  1.00 20.00      A    C  \n
HETATM   19  C26 KZN A   1       4.478   3.164   1.582  1.00 20.00      A    C  \n
HETATM   20  N27 KZN A   1       5.751   3.651  -0.394  1.00 20.00      A    N  \n
HETATM   21  O28 KZN A   1       6.622   3.940   0.095  1.00 20.00      A    O  \n
HETATM   22  N29 KZN A   1       0.447  -1.868   0.111  1.00 20.00      A    N  \n
HETATM   23  N30 KZN A   1       1.534  -1.125  -0.127  1.00 20.00      A    N  \n
HETATM   24  N32 KZN A   1      -3.791  -3.080  -1.856  1.00 20.00      A    N  \n
HETATM   25  C33 KZN A   1      -4.570  -3.213  -0.615  1.00 20.00      A    C  \n
HETATM   26  O34 KZN A   1      -5.679  -3.710  -0.645  1.00 20.00      A    O  \n
HETATM   27  O35 KZN A   1      -1.845  -2.378  -2.848  1.00 20.00      A    O  \n"""

	def __init__(self, pymolName):
		self.pymolName = pymolName


class TritylLabel:
	identifier = "T-R-I-T-Y-L"
	modifiedAA = True
	numberOfRotatingBonds = 8
	numberOfAtoms = 72
	rotate = True
	rotationInfo = {'1': [4,slice(6, numberOfAtoms), False], '2': [5,slice(7, numberOfAtoms), False], '3': [6,slice(8, numberOfAtoms), False], '4': [7,slice(9, numberOfAtoms), False], '5': [8,slice(10, numberOfAtoms), False], '6': [9,slice(11, numberOfAtoms), False], '7': [10,slice(12, numberOfAtoms), False], '8': [11,slice(13, numberOfAtoms), False] }
	clashAtoms = slice(5, numberOfAtoms)
	radius = 19.0
	atomNames =	 ['N', 'C', 'OAA', 'OAB' , 'CA', 'CB', 'SG', 'SD', 'CAD', 'CAE', 'OAC', 'CAF', 'CAG', 'CAH', 'CAI', 'CAJ', 'CAK','CAL', 'CAM', 'CAN', 'CAO', 'CAP', 'CAQ', 'CAR', 'CAS', 'CAT', 'CAU', 'CAV', 'CAW', 'CAX', 'CAY', 'CAZ', 'CBA', 'CBB', 'CBC', 'CBD', 'CBE', 'CBF', 'CBG', 'CBH', 'CBI', 'CBJ', 'CBK', 'CBL', 'CBM', 'CBN', 'CBO', 'CBP', 'CBQ', 'CBR', 'CBS', 'CBT', 'CBU', 'CBV', 'CBW', 'SAC', 'SAD', 'SAE', 'SAF', 'SAG', 'SAH', 'SAI', 'SAJ', 'SAK', 'SAL', 'SAM', 'SAN', 'OAD', 'OAE', 'OAF', 'OAG', 'OAH' ]
				#  0	   1	  2	     3	  4	   5	 6	   7	 8	   9	 10	   11	 12	  13	14	  15	 16	   17    18    19    20   21     22    23   24     25   26     27    28   29     30    31    32    33    34    35   36    37     38    39    40    41    42    43    44    45    46    47    48    49   50     51   52     53   54     55   56    57     58   59    60     61    62    63   64     65    66   67    68     69    70     71
	
	unsortedAtomNames = ['N', 'C', 'OAA', 'OAB' , 'CA', 'CB', 'SG', 'SD', 'CAD', 'CAE', 'OAC', 'CAF', 'CAG', 'CAH', 'CAI', 'CAJ', 'CAK','CAL', 'CAM', 'CAN', 'CAO', 'CAP', 'CAQ', 'CAR', 'CAS', 'CAT', 'CAU', 'CAV', 'CAW', 'CAX', 'CAY', 'CAZ', 'CBA', 'CBB', 'CBC', 'CBD', 'CBE', 'CBF', 'CBG', 'CBH', 'CBI', 'CBJ', 'CBK', 'CBL', 'CBM', 'CBN', 'CBO', 'CBP', 'CBQ', 'CBR', 'CBS', 'CBT', 'CBU', 'CBV', 'CBW', 'SAC', 'SAD', 'SAE', 'SAF', 'SAG', 'SAH', 'SAI', 'SAJ', 'SAK', 'SAL', 'SAM', 'SAN', 'OAD', 'OAE', 'OAF', 'OAG', 'OAH']
	movingAtoms=[]
	spinLocation = 'CAQ'
	highlight = 'CAQ'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "loose"
	numberToFind = {'painstaking': 2000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = "The calculations are very slow with this spin label!"
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 5.5
	pdbStr = """HEADER    TRITYL\n 	
ATOM      1  OAA   UNK   1       4.474   6.482  -4.055  1.00  0.00           O\n  
ATOM      2  C     UNK   1       5.495   6.703  -4.922  1.00  0.00           C\n  
ATOM      3  CA    UNK   1       6.756   5.901  -4.553  1.00  0.00           C\n  
ATOM      4  OAB   UNK   1       5.401   7.439  -5.888  1.00  0.00           O\n  
ATOM      5  N     UNK   1       7.938   6.334  -5.276  1.00  0.00           N\n  
ATOM      6  CB    UNK   1       6.461   4.399  -4.770  1.00  0.00           C\n  
ATOM      7  SG    UNK   1       7.735   3.261  -4.068  1.00  0.00           S\n  
ATOM      8  SD    UNK   1       7.668   3.655  -2.021  1.00  0.00           S\n  
ATOM      9  CAD   UNK   1       6.447   2.445  -1.349  1.00  0.00           C\n  
ATOM     10  CAE   UNK   1       5.029   2.692  -1.839  1.00  0.00           C\n  
ATOM     11  OAC   UNK   1       4.155   1.779  -1.141  1.00  0.00           O\n  
ATOM     12  CAF   UNK   1       2.867   1.743  -1.576  1.00  0.00           C\n  
ATOM     13  CAG   UNK   1       2.029   0.763  -0.846  1.00  0.00           C\n  
ATOM     14  OAD   UNK   1       2.454   2.477  -2.466  1.00  0.00           O\n  
ATOM     15  CAH   UNK   1       0.616   0.906  -0.952  1.00  0.00           C\n  
ATOM     16  CAI   UNK   1       2.556  -0.304  -0.062  1.00  0.00           C\n  
ATOM     17  SAC   UNK   1      -0.151   2.197  -1.896  1.00  0.00           S\n  
ATOM     18  CAJ   UNK   1      -0.252  -0.011  -0.311  1.00  0.00           C\n  
ATOM     19  SAD   UNK   1       4.294  -0.625   0.148  1.00  0.00           S\n  
ATOM     20  CAK   UNK   1       1.685  -1.190   0.613  1.00  0.00           C\n  
ATOM     21  CAL   UNK   1      -1.808   2.056  -1.045  1.00  0.00           C\n  
ATOM     22  SAE   UNK   1      -1.987   0.275  -0.544  1.00  0.00           S\n  
ATOM     23  CAM   UNK   1       0.261  -1.086   0.485  1.00  0.00           C\n  
ATOM     24  CAN   UNK   1       4.074  -2.399   0.674  1.00  0.00           C\n  
ATOM     25  SAF   UNK   1       2.464  -2.442   1.603  1.00  0.00           S\n  
ATOM     26  CAO   UNK   1      -2.909   2.399  -2.056  1.00  0.00           C\n  
ATOM     27  CAP   UNK   1      -1.842   2.965   0.189  1.00  0.00           C\n  
ATOM     28  CAQ   UNK   1      -0.633  -2.045   1.142  1.00  0.00           C\n  
ATOM     29  CAR   UNK   1       5.212  -2.782   1.625  1.00  0.00           C\n  
ATOM     30  CAS   UNK   1       4.003  -3.314  -0.554  1.00  0.00           C\n  
ATOM     31  CAT   UNK   1      -1.813  -2.555   0.435  1.00  0.00           C\n  
ATOM     32  CAU   UNK   1      -0.348  -2.497   2.515  1.00  0.00           C\n  
ATOM     33  CAV   UNK   1      -1.755  -2.916  -0.950  1.00  0.00           C\n  
ATOM     34  CAW   UNK   1      -3.072  -2.704   1.103  1.00  0.00           C\n  
ATOM     35  CAX   UNK   1      -0.045  -1.566   3.559  1.00  0.00           C\n  
ATOM     36  CAY   UNK   1      -0.364  -3.884   2.857  1.00  0.00           C\n  
ATOM     37  SAG   UNK   1      -0.267  -2.819  -1.912  1.00  0.00           S\n  
ATOM     38  CAZ   UNK   1      -2.902  -3.386  -1.634  1.00  0.00           C\n  
ATOM     39  SAH   UNK   1      -3.285  -2.352   2.832  1.00  0.00           S\n  
ATOM     40  CBA   UNK   1      -4.232  -3.124   0.411  1.00  0.00           C\n  
ATOM     41  SAI   UNK   1      -0.201   0.184   3.323  1.00  0.00           S\n  
ATOM     42  CBB   UNK   1       0.240  -1.994   4.878  1.00  0.00           C\n  
ATOM     43  SAJ   UNK   1      -0.493  -5.154   1.632  1.00  0.00           S\n  
ATOM     44  CBC   UNK   1      -0.106  -4.322   4.179  1.00  0.00           C\n  
ATOM     45  CBD   UNK   1      -0.833  -4.005  -3.227  1.00  0.00           C\n  
ATOM     46  SAK   UNK   1      -2.688  -3.813  -3.342  1.00  0.00           S\n  
ATOM     47  CBE   UNK   1      -4.155  -3.482  -0.967  1.00  0.00           C\n  
ATOM     48  CBF   UNK   1      -5.108  -2.021   2.701  1.00  0.00           C\n  
ATOM     49  SAL   UNK   1      -5.729  -3.168   1.370  1.00  0.00           S\n  
ATOM     50  CBG   UNK   1       0.858   0.632   4.781  1.00  0.00           C\n  
ATOM     51  SAM   UNK   1       0.551  -0.699   6.054  1.00  0.00           S\n  
ATOM     52  CBH   UNK   1       0.202  -3.382   5.203  1.00  0.00           C\n  
ATOM     53  CBI   UNK   1      -0.975  -6.485   2.835  1.00  0.00           C\n  
ATOM     54  SAN   UNK   1      -0.103  -6.074   4.442  1.00  0.00           S\n  
ATOM     55  CBJ   UNK   1      -0.213  -3.595  -4.567  1.00  0.00           C\n  
ATOM     56  CBK   UNK   1      -0.479  -5.445  -2.836  1.00  0.00           C\n  
ATOM     57  CBL   UNK   1      -5.320  -3.938  -1.762  1.00  0.00           C\n  
ATOM     58  CBM   UNK   1      -5.781  -2.389   4.028  1.00  0.00           C\n  
ATOM     59  CBN   UNK   1      -5.355  -0.561   2.306  1.00  0.00           C\n  
ATOM     60  CBO   UNK   1       0.386   1.974   5.351  1.00  0.00           C\n  
ATOM     61  CBP   UNK   1       2.339   0.651   4.385  1.00  0.00           C\n  
ATOM     62  CBQ   UNK   1       0.493  -3.918   6.555  1.00  0.00           C\n  
ATOM     63  CBR   UNK   1      -0.453  -7.834   2.326  1.00  0.00           C\n  
ATOM     64  CBS   UNK   1      -2.493  -6.485   3.047  1.00  0.00           C\n  
ATOM     65  OAE   UNK   1      -6.361  -4.342  -0.989  1.00  0.00           O\n  
ATOM     66  OAF   UNK   1      -5.341  -3.976  -2.988  1.00  0.00           O\n  
ATOM     67  OAG   UNK   1       0.489  -2.960   7.515  1.00  0.00           O\n  
ATOM     68  OAH   UNK   1       0.707  -5.103   6.791  1.00  0.00           O\n  
ATOM     69  CBT   UNK   1      -7.578  -4.731  -1.671  1.00  0.00           C\n  
ATOM     70  CBU   UNK   1       0.818  -3.367   8.866  1.00  0.00           C\n 
ATOM     71  CBV   UNK   1      -8.455  -3.526  -1.987  1.00  0.00           C\n  
ATOM     72  CBW   UNK   1      -0.428  -3.764   9.645  1.00  0.00           C\n"""

	
	def __init__(self, pymolName):
		self.pymolName = pymolName

class RxLabel:
	identifier = "R-X"
	modifiedAA = True
	numberOfRotatingBonds = 8
	numberOfAtoms = 72
	rotate = True
	rotationInfo = {}
	clashAtoms = slice(5, numberOfAtoms)
	radius = 7.5
	atomNames = ['N1']
	unsortedAtomNames = ['N1']
	movingAtoms=[]
	spinLocation = 'N1'
	highlight = 'N1'
	atomsForSuperposition = []
	defaultVdw = "loose"
	#numberToFind = {'painstaking': 2000, 'thorough search': 200, 'quick search': 50}
	#numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = "The location of this label is calculated by a simpler algorithm. No rotamers are created!"
	errorMessage = ""
	trialAtomSphereRadius = 7.5
	exclusionSphereRadius = 6.5
	pdbStr = """HEADER    RX\n 	
ATOM      1 N1   RXR A   1       8.144   2.482  -5.989  1.00 20.00           N\n\n"""

	
	def __init__(self, pymolName):
		self.pymolName = pymolName
