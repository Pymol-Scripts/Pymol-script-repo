"""
---mtsslWizard: spin_labeling plugin for PyMOL --- 
Author	: Gregor Hagelueken
Date	: Jan 2013
Version : 1.1
Mail	: hagelueken'at'pc.uni-bonn.de
 
mtsslWizard is a plugin for the PyMOL Molecular Graphics System. 
It allows in silico spin labeling of proteins with the spin label MTSSL in PyMOL. Also, distances between ensembles of two spin labels can be calculated and exported.
The program was tested with PyMOL version 1.5.
 
Please cite:
Hagelueken G, Ward R, Naismith JH, Schiemann O. MtsslWizard: In silico Spin-Labeling and Generation of Distance Distributions in PyMOL. 2012. Appl. Mag. Res., accepted for publication.
 
----------------------------------------------------------------------
----------------------------------------------------------------------
 
"""
import pymol
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
from Tkinter import Tk

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

##########################
#wizard class            #
##########################
class MtsslWizard(Wizard):
	def __init__(self):
		#print platform.system()
		Wizard.__init__(self)
		print "**************************************************************************************************"
		print "* MtsslWizard by gha.                                                                            *"
		print "* Please remove any solvent or unwanted heteroatoms before using the wizard!                     *"
		print "* You can do this e.g. by issuing 'remove solvent'.                                              *"
		print "**************************************************************************************************"
		
		#create contents of wizard menu
		self.reset()
		self.menu['mode'] = [
									  [1, 'Search','cmd.get_wizard().set_mode("Search")'],
									  [1, 'Measure','cmd.get_wizard().set_mode("Measure")'],
									  [1, 'Copy & Move','cmd.get_wizard().set_mode("Copy & Move")']
									  ]
		self.menu['currentLabel'] = [
									  [ 2, '\\559Protein', ''],
									  [1, 'MTSSL','cmd.get_wizard().set_currentLabel("MTSSL")'],
									  [1, 'PROXYL','cmd.get_wizard().set_currentLabel("PROXYL")'],
									  [ 2, '\\559DNA', ''],
									  [1, 'URIP','cmd.get_wizard().set_currentLabel("URIP")'],
									  [1, 'C','cmd.get_wizard().set_currentLabel("CLABEL")'],
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
 
		self.menu['clashes'] = [
									  [ 2, '\\955NOTE:\\559This parameter should normally ', ''],
									  [ 2, '\\559be set to "0"! Use it with care and only ', ''],
									  [ 2, '\\559change it if conformational changes of ', ''],
									  [ 2, '\\559the protein are expected!!!', ''],
									  [1, '0 clashes','cmd.get_wizard().set_clashes(0)'],
									  [1, '1 clashes','cmd.get_wizard().set_clashes(1)'],
									  [1, '2 clashes','cmd.get_wizard().set_clashes(2)'],
									  [1, '3 clashes','cmd.get_wizard().set_clashes(3)'],
									  [1, '4 clashes','cmd.get_wizard().set_clashes(4)'],
									  [1, '5 clashes','cmd.get_wizard().set_clashes(5)']
									  ]
 
		self.menu['cutoff'] = [
									  [ 2, '\\955NOTE:\\559This parameter determines ', ''],
									  [ 2, '\\559the minimum allowed distance between', ''],
									  [ 2, '\\559label and	protein atoms. Only change', ''],
									  [ 2, '\\559this parameter if conformational ch-', ''],
									  [ 2, '\\559anges at the labeling site are expected!', ''],
									  [1, '2.6 A','cmd.get_wizard().set_cutoff(2.6)'],
									  [1, '2.8 A','cmd.get_wizard().set_cutoff(2.8)'],
									  [1, '3.0 A','cmd.get_wizard().set_cutoff(3.0)'],
									  [1, '3.2 A','cmd.get_wizard().set_cutoff(3.2)'],
									  [1, '3.4 A','cmd.get_wizard().set_cutoff(3.4)'],
									  [1, '3.6 A','cmd.get_wizard().set_cutoff(3.6)'],
									  [1, '3.8 A','cmd.get_wizard().set_cutoff(3.8)']
									  ]
		
		self.menu['writeToFile'] = [
									  [1, 'no','cmd.get_wizard().set_writeToFile("no")'],
									  [1, 'yes','cmd.get_wizard().set_writeToFile("yes")']
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
		def urip():
			tmp=UripLabel("tmp")
			return tmp
		def clabel():
			tmp=CLabel("tmp")
			return tmp
		options = {"MTSSL":mtssl, "PROXYL":proxyl, "URIP":urip, "CLABEL":clabel,}
		self.set_vdwRestraints(options[currentLabel]().defaultVdw)
		print options[currentLabel]().info
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
		if self.mode == 'Search' or 'Stochastic':
			self.prompt = [ 'Select a residue to label...']
		if self.pick_count == 0 and self.mode == 'Measure':
			self.prompt = [ 'Select first label...']
		if self.pick_count == 0 and self.mode == 'Copy & Move':
			self.prompt = [ 'Select label to be copied & moved ...']
		if self.pick_count == 1 and self.mode == 'Measure':
			self.prompt = [ 'Select second label...']
		if self.pick_count == 1 and self.mode == 'Copy & Move':
			self.prompt = [ 'Select position for copied label ...']
		if self.running:
			self.prompt = [ 'Running, please wait...' ]
		return self.prompt
	
	def get_panel(self):
		if not ((self.mode == 'Measure') or (self.mode == 'Copy & Move')):
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
		elif self.mode == 'Copy & Move':
			return [
					[ 1, 'Mtssl Wizard',''],
					[ 3, 'Mode: %s'%self.mode,'mode'],
					[ 2, 'Reset','cmd.get_wizard().reset()'],
					[ 2, 'Done','cmd.set_wizard()']
					]
	#reset wizard to defaults
	def reset(self):
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
		self.writeToFile='no'
		self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
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
		print "Deleting everything..."
		cmd.delete(self.object_prefix + "*")
 
	def delete_last(self):
		try:
			print self.numberOfLabel
			if self.numberOfLabel >= 1:
				cmd.delete(self.object_prefix+str(self.numberOfLabel)+"*")
				self.numberOfLabel-=1
		except pymol.CmdException, pmce:
			print pmce
	
	def cleanup(self):
		print "Cleaning up..."
		cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
		#self.reset()
		#self.delete_all()
 
	def do_select(self, name):
		try:
			self.do_pick(0)
		except pymol.CmdException, pmce:
			print pmce
 
	def do_pick(self, picked_bond):
		#This is to catch repeated selections in Search mode
		if self.pick_count > 0 and self.mode == 'Search':
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
				print "MtsslWizard: object not found."
			self.pick_count += 1
			if self.mode == 'Measure' or self.mode == 'Copy & Move':
				#deselect before next pick
				cmd.deselect()
			self.cmd.refresh_wizard()

		#second click	 
		elif self.pick_count == 1  and (self.mode == 'Measure' or self.mode == 'Copy & Move'):
			self.residue2_name = self.createSelectionMacro("(sele)")
			#print self.residue2_name
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
				print "MtsslWizard: object not found."
			self.pick_count += 1
			#deselect before next pick
			self.run()
			cmd.deselect()
		
		if self.mode == 'Search' or self.mode == 'Stochastic':
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
		
	##########################
	#RUN                     #
	##########################
	def run(self):
		self.conformationList = []
		my_view = cmd.get_view()
		##########################
		#Search mode             #
		##########################
		if self.pick_count == 1 and self.mode == 'Search':
			print "\n\n\nNew run:\n"
			#generate the label and superpose onto selected position
			if self.currentLabel == "MTSSL":
				self.label = MtsslLabel("mtssl_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "PROXYL":
				self.label = ProxylLabel("proxyl_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "URIP":
				self.label = UripLabel("urip_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			elif self.currentLabel == "CLABEL":
				self.label = CLabel("C_"+str(self.numberOfLabel))
				cmd.read_pdbstr(self.label.pdbStr, self.label.pymolName)
			print self.residue1_name, self.label.identifier
			
			#attach identifier for identification in pymol
			#cmd.set_name(self.residue1_name, "%s_%s" %(self.residue1_name, self.label.identifier))
			#self.residue1_name = "%s_%s" %(self.residue1_name, self.label.identifier)
			
			print "Attempting superposition..."
			if not self.superpose():
				print "Superposition does not work."
				print "Possible reasons:"
				print "1) Glycine? Mutate to Ala first."
				print "2) Trying to attach DNA label to Protein or vice versa?"
				if len(self.label.errorMessage) > 0:
					print "3) %s" %self.label.errorMessage
				self.cleanupAfterRun(my_view)
				return
			else:
				print "Superposition worked!"
			
			#prepare movingAtoms, put into correct order...
			stored.movingAtoms = []
			cmd.iterate_state(1, self.label.pymolName, 'stored.movingAtoms.append((x,y,z))')
			atoms = stored.movingAtoms
			self.label.prepareMovingAtoms(atoms)
			
			#create object with only the atoms around the label to speed everything up 
			#cmd.color ("red", "%s &! %s within %f of %s" %(self.picked_object1, self.residue1_name, self.label.radius, self.label.pymolName))
			protein="%s &! %s within %f of %s" %(self.picked_object1, self.residue1_name, self.label.radius, self.label.pymolName)
			stored.environmentAtoms = []
			cmd.iterate_state(1, protein, 'stored.environmentAtoms.append((x,y,z))')
			environmentAtoms=numpy.array(stored.environmentAtoms)
			result = self.fastMtsslify(environmentAtoms)
			#only switch on snuggly fit search for "painstaking"
			if self.thoroughness == "painstaking":
				self.createSnugglyFitConformations()		
			print ""
			print "Found: %i in %i tries." %(result[0], result[1])
			if result[0] > 0 and result[0] <= 10 and self.vdwRestraints == "tight" and not self.currentLabel == "CLABEL":
				print "The number of conformations is very small!"
				print "Consider switching vdW restraints to 'loose' for this position!"
			print "Done!"
			self.finalCosmetics(result[0])
			
		##########################
		#Measure mode            #
		##########################		
		elif self.pick_count == 2 and self.mode == 'Measure':
			print "\n\n\nDistance calculation:\n"
			print "The dashed lines are the c-beta distance (green),\nand the distance between the geometric averages\nof the two ensembles (yellow).\n"
			print "The following statistics refer to the distribution\nof the individual distances between all conformers (may take a while):\n"
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
			def urip():
				tmp=UripLabel("tmp")
				return tmp.spinLocation
			def clabel():
				tmp=CLabel("tmp")
				return tmp.spinLocation
			options = {"M-T-S-S-L":mtssl, "P-R-O-X-Y-L":proxyl, "U-R-I-P":urip, "C-L-A-B-E-L":clabel,}
			
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
			cmd.distance(self.object_prefix+"avg","tmp_average1 & name ps1","tmp_average2 & name ps1")
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
			
			#Copy to clipboard
			self.copyStringToClipboard(outputStr)

			#Write to file
			if self.writeToFile=='yes':
				try:
					filename = "%s-%s" %(self.residue1_name, self.residue2_name)
					numpy.savetxt(filename, output, delimiter='\t')
					print "Written to file:"
					print "%s/%s" %(os.getcwd(), filename)
				except:
					print "Writing to file failed!"
					
			print calculateStatistics2(dist)
			if len(cBeta) > 0:
				print "Cbeta distance: %3.1f" %cBeta[0]
		
		##########################
		#copy & move mode        #
		##########################	
		elif self.pick_count == 2 and self.mode == 'Copy & Move':
			print "\n\n\nCopy labels:\n"
			self.copyAndMove()
				
		#some cleanup after each run
		cmd.set_view(my_view)
		self.cleanupAfterRun(my_view)
	
	##########################
	#various methods         #
	##########################
	def copyStringToClipboard(self, string):
		try:
			r = pymol._ext_gui.root
			r.clipboard_clear()
			r.clipboard_append(string)
			print "Copied to clipboard."
			return
		except:
			pass
		try:
			import pyperclip
			pyperclip.copy(string)
			print "Copied to clipboard."
			return
		except:
			pass
		try:
			import xerox
			xerox.copy(string)
			print "Copied to clipboard."
			return
		except:
			pass
		print "Copy to clipboard failed. Try to install either the 'pyperclip' or 'xerox' module for Python."
	
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
	
	##########################
	#superpose               #
	##########################	
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
		if apply(cmd.pair_fit, args):
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
	
	##########################
	#fastMtsslify            #
	##########################	
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
		print "Trying to find conformations for label %s with vdW restraints: %s" % (self.label.pymolName, self.vdwRestraints)
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
					print found,
					#only switch on snuggly fit search for "painstaking"
					if self.thoroughness == "painstaking":
						vdw = numberOfVdwContacts(self.label.movingAtoms[self.label.clashAtoms], environmentAtoms, self.cutoff)
						thisConformation = [found, vdw]
						self.conformationList.append(thisConformation)	
					self.createRotamer(found)
				else:
					print "i",
			else:
				print ".",
			ntries+=1
		results = [found, ntries]
		return results
	
	##########################
	#createRotamer           #
	##########################
	def createRotamer(self, found):
		for i in range (0, len(self.label.movingAtoms)):
			stored.xyz = []
			stored.xyz = self.label.movingAtoms[i]
			cmd.alter_state(1,self.label.pymolName +"& name " + self.label.atomNames[i],"(x,y,z)=stored.xyz")
		cmd.create("%s_%s" %(self.residue1_name, self.label.identifier), self.label.pymolName, 1, found)
	
	################################
	#createSnugglyFitConformations #
	################################
	def createSnugglyFitConformations(self): #select conformations with the best fit to the molecular surface, rank them and create an object
		print "Snuggliest fit(s):"
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
				print "Conformation %i: \t\t%i \t\t\t vdW contacts" %(snugglyFitList[i]['state'],snugglyFitList[i]['vdw'])
		print "Average vdW contacts of all possible conformations: ",averageAtomCount
	
	##########################
	#final cosmetics         #
	##########################
	def finalCosmetics(self, found): #make everything look nice
		if found >= 1:
			print "Found %i conformations." %found
			#show all conformations and do some coloring
			cmd.set("all_states",1)
			self.toggleStatesCaption='Toggle states: ON'
			cmd.color("blue","%s_%s" %(self.residue1_name, self.label.identifier))
			cmd.color("red", "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.highlight))
			util.cnc("%s_%s" %(self.residue1_name, self.label.identifier))
			cmd.disable(self.label.pymolName) 
			cmd.show("spheres", "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.highlight))
			cmd.set("sphere_scale", "0.2", "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.highlight))
			identifierLabel="%s|%s|%s|%s" %(self.residue1_name, self.currentLabel, self.vdwRestraints, self.thoroughness)
			print identifierLabel
			#mark average N1 position with pseudoatom
			stored.label = []
			cmd.iterate_state(0, "%s_%s & name %s" %(self.residue1_name, self.label.identifier, self.label.spinLocation), 'stored.label.append((x,y,z))')
			atoms1=numpy.array(stored.label)
			#create pseudoatom at average coordinate of each ensemble
			avgAtoms=numpy.average(atoms1,axis=0)
			self.createPseudoatom (avgAtoms, "%s_label" %(self.residue1_name), 1)
			cmd.set("sphere_scale", "0.5", "%s_label" %(self.residue1_name))
			cmd.label("%s_label" %(self.residue1_name), `identifierLabel`)
			cmd.show("label")
			cmd.show("spheres", "name PS1")
			#print label+"*"+","+"labelEnvironment_"+label
			cmd.delete(self.label.pymolName+"*")
			cmd.group("%s%s" %(self.object_prefix, str(self.numberOfLabel)), "%s*, labelEnvironment_%s,%s*" %(self.label.pymolName, self.label.pymolName, self.residue1_name))
		else:
			print "Sorry, I did not find anything. Your options are:\n1) Try again, maybe with increased thoroughness,\n2) change vdW restraints to 'loose'"
			cmd.delete("%s*" %(self.label.pymolName))

##########################
#start wizard            #
##########################						  
def open_wizard():
	wiz = MtsslWizard()
	cmd.set_wizard(wiz)

##########################
#helper functions        #
##########################
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
		duplicated = True
		atoms1=numpy.tile(atoms1, (2,1))
	if len(atoms2) == 1:
		duplicated = True
		atoms2=numpy.tile(atoms2, (2,1))
		
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
	spinLocation = 'N1'
	highlight = 'O1'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 10000, 'quick search': 1000}
	info = ""
	errorMessage = ""
	pdbStr = """HEADER    MTSSL\n
COMPND    coordinates of R1A      from program: libcheck\n                        
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1\n                      
SCALE1      0.010000  0.000000  0.000000        0.00000\n                         
SCALE2      0.000000  0.010000  0.000000        0.00000\n                         
SCALE3      0.000000  0.000000  0.010000        0.00000\n                         
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
	
		
	def prepareMovingAtoms(self, atoms):
		movingAtoms = atoms
		#This is the order of atoms when iterate is used on mtssl in pymol:
		#['N', 'CA', 'C', 'O', 'CB', 'SG', 'SD', 'CE', 'N1', 'O1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#  0	1	 2	  3		4	  5		6	  7		8	  9		10	  11	12	  13	14	  15	16	  17  
		#This is the order we want
		#['N', 'O', 'C', 'CA', 'CB', 'SG', 'SD', 'CE', 'C3', 'O1', 'C2', 'N1', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#  0	1	 2	  3		4	  5		6	  7		8	  9		10	  11	12	  13	14	  15	16	  17  
		
		#swap N1 with C3 and CA with O, so that all atoms which serve as rotation axes are in sequence
		n1=movingAtoms[8]
		c3=movingAtoms[11]
		ca=movingAtoms[1]
		o=movingAtoms[3]
		movingAtoms[8]=c3
		movingAtoms[11]=n1
		movingAtoms[1]=o
		movingAtoms[3]=ca
		self.movingAtoms=numpy.array(movingAtoms)

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
	unsortedAtomNames = ['N', 'CA', 'C', 'O', 'CB', 'SG', 'C1', 'N1', 'O1', 'C2', 'N2', 'C3', 'O3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10']
	spinLocation = 'N1'
	highlight = 'O1'
	atomsForSuperposition = ['CA','N','C','CB']
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 20000, 'quick search': 3000}
	info = ""
	errorMessage = ""
	pdbStr = """HEADER    PROXYL\n
COMPND    coordinates of PROXYL   from program: MMM\n
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1\n                      
SCALE1      0.010000  0.000000  0.000000        0.00000\n                         
SCALE2      0.000000  0.010000  0.000000        0.00000\n                         
SCALE3      0.000000  0.000000  0.010000        0.00000\n                         
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
	
		
	def prepareMovingAtoms(self, atoms):
		movingAtoms = atoms
		#This is the order of atoms when iterate is used on proxyl in pymol:
		#['N', 'CA', 'C', 'O', 'CB', 'SG', 'C1', 'N1', 'O1', 'C2', 'N2', 'C3', 'O3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10']
		#  0	1	  2	   3	4	  5		6	  7		8	  9		10	  11	12	  13	14	  15	16	  17	 18	   19 
		#This is the order we want:
		#['N', 'O', 'C', 'CA', 'CB', 'SG', 'C1', 'C2', 'N2', 'C3', 'N1', 'O1', 'C10','C4', 'C5', 'C6', 'C7', 'C8', 'C9',  'O3']
		#  0	1	 2	  3		4	  5		6	  7		8	  9		10	  11	12	  13	14	  15	16	  17	 18	   19 
		
		#swap atoms
		ca=movingAtoms[1]
		o=movingAtoms[3]
		n1=movingAtoms[7]
		o1=movingAtoms[8]
		c2=movingAtoms[9]
		n2=movingAtoms[10]
		c3=movingAtoms[11]
		o3=movingAtoms[12]
		c10=movingAtoms[19]
		
		movingAtoms[1]=o
		movingAtoms[3]=ca
		movingAtoms[7]=c2
		movingAtoms[8]=n2
		movingAtoms[9]=c3
		movingAtoms[10]=n1
		movingAtoms[11]=o1
		movingAtoms[12]=c10
		movingAtoms[19]=o3
		self.movingAtoms=numpy.array(movingAtoms)
		
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
	unsortedAtomNames = ["O4'", "C3'", "C2'", 'CS1', 'CS10', 'CS11', 'CS2', 'CS3', 'CS5', 'CS6', 'CS7', 'CS8', 'CS9', 'NS1', 'NS2', 'NS3', 'OS1', 'OS2']
	spinLocation = 'NS3'
	highlight = 'OS1'
	atomsForSuperposition = ["C2'","C3'","O4'"]
	defaultVdw = "loose"
	numberToFind = {'painstaking': 1000, 'thorough search': 200, 'quick search': 50}
	numberOfTries = {'painstaking': 100000, 'thorough search': 20000, 'quick search': 3000}
	info = "\nFor the Urip label, the vdW restraints are by default set to 'loose' to account for possible polar interactions between the amide bonds and the DNA backbone."
	errorMessage = "Check atom nomenclature. The ribose atoms are sometimes called C2* instead of C2'\nThis can be changed e.g. by 'alter all, name=string.replace(name,\"*\",\"'\")' in PyMOL."
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
	
		
	def prepareMovingAtoms(self, atoms):
		movingAtoms = atoms
		#This is the order of atoms when iterate is used on urip in pymol:
		#["O4'", "C3'", "C2'", 'CS1', 'CS10', 'CS11', 'CS2', 'CS3', 'CS5', 'CS6', 'CS7', 'CS8', 'CS9', 'NS1', 'NS2', 'NS3', 'OS1', 'OS2']
		#  0	   1	 2	    3	   4	   5	   6	  7		 8	    9	   10	  11	 12	    13	   14	  15	 16	    17 
		#This is the order we want:
		#["O4'", "C3'", "C2'", 'NS1', 'CS1', 'NS2', 'CS7', 'CS2', 'CS10', 'CS11',  'CS3', 'CS5', 'CS6',  'CS8', 'CS9',   'NS3', 'OS1', 'OS2']
		#  0	  1	     2	    3	   4	  5		 6	    7	   8	   9		10	   11	  12	  13	 14	      15	 16	    17
		#swap atoms
		cs1=movingAtoms[3]
		cs10=movingAtoms[4]
		cs11=movingAtoms[5]
		cs2=movingAtoms[6]
		cs3=movingAtoms[7]
		cs5=movingAtoms[8]
		cs6=movingAtoms[9]
		cs7=movingAtoms[10]
		cs8=movingAtoms[11]
		cs9=movingAtoms[12]
		ns1=movingAtoms[13]
		ns2=movingAtoms[14]
		
		movingAtoms[3]=ns1
		movingAtoms[4]=cs1
		movingAtoms[5]=ns2
		movingAtoms[6]=cs7
		movingAtoms[7]=cs2
		movingAtoms[8]=cs10
		movingAtoms[9]=cs11
		movingAtoms[10]=cs3
		movingAtoms[11]=cs5
		movingAtoms[12]=cs6
		movingAtoms[13]=cs8
		movingAtoms[14]=cs9
		
		
		
		self.movingAtoms=numpy.array(movingAtoms)

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
	unsortedAtomNames = ['C1', 'N1', 'C2', 'O2', 'N3', 'C5', 'C6', 'N7', 'C8', 'C9', 'O10', 'C11', 'C12', 'C13', 'C14', 'C15', 'N16', 'C17', 'C19', 'C20', 'O21', 'C22', 'C23']
	spinLocation = 'N16'
	highlight = 'O21'
	atomsForSuperposition = ["N1", "C2", "O2", "N3"]
	#atomsForSuperposition = ["C2'","C3'","O4'"]
	defaultVdw = "tight"
	numberToFind = {'painstaking': 1, 'thorough search': 1, 'quick search': 1}
	numberOfTries = {'painstaking': 1, 'thorough search': 1, 'quick search': 1}
	info = "\nThe C label is only superimposed onto e.g. dC but does not move.\n"
	errorMessage = "This label does not superpose onto A or G!"
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
	
		
	def prepareMovingAtoms(self, atoms):
		movingAtoms = atoms
		self.movingAtoms=numpy.array(movingAtoms)