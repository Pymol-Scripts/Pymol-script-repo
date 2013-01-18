"""
---mtsslWizard: spin_labeling plugin for PyMOL --- 
Author  : Gregor Hagelueken
Date    : Jan 2013
Version : 1.1
Mail    : hagelueken'at'pc.uni-bonn.de
 
mtsslWizard is a plugin for the PyMOL Molecular Graphics System. 
It allows in silico spin labeling of proteins with the spin label MTSSL in PyMOL. Also, distances between ensembles of two spin labels can be calculated and exported.
The program was tested with PyMOL versions 1.4.
 
Please cite:
 Hagelueken G, Ward R, Naismith JH, Schiemann O. MtsslWizard: In silico Spin-Labeling and Generation of Distance Distributions in PyMOL. 2012. Appl. Mag. Res., accepted for publication.
 
----------------------------------------------------------------------
----------------------------------------------------------------------
 
"""
import pymol
import numpy
import scipy.spatial.distance
import random, time, math
from pymol import cmd
from pymol.wizard import Wizard
from pymol import stored
from operator import itemgetter
from Tkinter import Tk

default_thoroughness = "normal search"
default_label = "MTSSL"
default_cutoff = 3.4
default_clashes = 0
default_mode = 'Search'
default_rotamers = 'Unrestricted'
default_clashguard = 'on'
internalClash_cutoff = 2.5

def __init__(self):
    self.menuBar.addmenuitem('Wizard', 'command',
                             'MtsslWizard',
                             label = 'MtsslWizard',
                             command = lambda s=self : open_wizard())
class MtsslWizard(Wizard):
    def __init__(self):
        Wizard.__init__(self)
        print "MtsslWizard by gha. Please remove any solvent or unwanted heteroatoms before using the wizard!"
        print "!!! This is an alpha version for testing purposes. Please do not distribute it!!!"
        
        #create contents of wizard menu
        self.reset()
        self.menu['mode'] = [
                                      [1, 'Search','cmd.get_wizard().set_mode("Search")'],
                                      [1, 'Measure','cmd.get_wizard().set_mode("Measure")'],
                                      [1, 'Copy & Move','cmd.get_wizard().set_mode("Copy & Move")']
                                      ]
        self.menu['currentLabel'] = [
                                      [1, 'MTSSL','cmd.get_wizard().set_currentLabel("MTSSL")'],
                                      [1, 'PROXYL','cmd.get_wizard().set_currentLabel("PROXYL")'],
                                      ]
        self.menu['thoroughness'] = [
                                      [1, 'painstaking','cmd.get_wizard().set_thoroughness("painstaking")'],
                                      [1, 'thorough search','cmd.get_wizard().set_thoroughness("thorough search")'],
                                      [1, 'normal search','cmd.get_wizard().set_thoroughness("normal search")'],
                                      [1, 'quick search','cmd.get_wizard().set_thoroughness("quick search")']
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
                                      [ 2, '\\559label and  protein atoms. Only change', ''],
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

    #some setter functions    
    def set_rotamers(self,rotamers):
        self.rotamers = rotamers
        self.cmd.refresh_wizard()

    def set_currentLabel(self, currentLabel):
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
    
    #reset wizard to defaults
    def reset(self):
        self.toggleStatesCaption='Toggle states: ON'
        self.start_time=time.time()
        self.conformationList=[]
        self.object_prefix = "mW-"
        self.pick_count = 0
        self.object_count = 0
        self.allowedAngle=[False,False,False,False,False]
        self.thoroughness = default_thoroughness
        self.cutoff = default_cutoff
        self.clashes = default_clashes
        self.rotamers = default_rotamers
        self.currentLabel = default_label
        self.residue1_name = None
        self.residue2_name = None
        self.picked_object1 = None
        self.picked_object2 = None
        self.numberOfLabel = 0
        self.mode=default_mode
        self.writeToFile='no'
        self.clashGuard=default_clashguard
        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
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
        
    def set_mode(self, mode):
        self.mode = mode
        self.cmd.refresh_wizard()
 
    def cleanup(self):
        print "Cleaning up..."
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        #self.reset()
        #self.delete_all()
 
    def get_prompt(self):
        self.prompt = None
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
        return self.prompt
 
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
            self.residue1_name = createSelectionMacro("(sele)") 
            # transfer the click selection to a named selection
            cmd.select(self.residue1_name+"_tmp", "(sele)")
            # find the name of the object which contains the selection
            new_name = None
            obj_list = cmd.get_names('objects')
            for object in obj_list:
                if cmd.get_type(object)=="object:molecule":
                    if cmd.count_atoms("(%s and %s)"%(object, self.residue1_name+"_tmp")):
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
            self.residue2_name = createSelectionMacro("(sele)")
            #print self.residue2_name
            # transfer the click selection to a named selection
            cmd.select(self.residue2_name+"_tmp", "(sele)")

            # find the name of the object which contains the selection
            new_name = None
            obj_list = cmd.get_names('objects')
            for object in obj_list:
                if cmd.get_type(object)=="object:molecule":
                    if cmd.count_atoms("(%s and %s)"%(object, self.residue2_name+"_tmp")):
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
		
    def run(self):
        self.conformationList=[]
        my_view= cmd.get_view()
        
        #Search mode
        if self.pick_count == 1 and self.mode == 'Search':
            print "\n\n\nNew run:\n"
            fastMtsslify(self, self.picked_object1, self.residue1_name, self.numberOfLabel, self.thoroughness, self.cutoff, self.clashes)
        
        #Measure mode       
        elif self.pick_count == 2 and self.mode == 'Measure':
			print "\n\n\nDistance calculation:\n"
			print "The dashed lines are the c-beta distance (green),\nand the distance between the geometric averages\nof the two ensembles (yellow)."
			print "The following statistics refer to the distribution\nof the individual distances between all conformers:\n"
			stored.label1 = []
			stored.label2 = []
			
			#account for measuring to single atoms, such as metal centers and different labels
			numberOfAtomsInS1=cmd.count_atoms(self.residue1_name)
			numberOfAtomsInS2=cmd.count_atoms(self.residue2_name)
			if numberOfAtomsInS1==18 or numberOfAtomsInS1==20:
				cmd.iterate_state(0, self.residue1_name+" & name N1", 'stored.label1.append((x,y,z))')
			else:
				cmd.iterate_state(0, self.residue1_name, 'stored.label1.append((x,y,z))')
			if numberOfAtomsInS2==18  or numberOfAtomsInS1==20:
				cmd.iterate_state(0, self.residue2_name+" & name N1", 'stored.label2.append((x,y,z))')
			else:
				cmd.iterate_state(0, self.residue2_name, 'stored.label2.append((x,y,z))')

			atoms1=numpy.array(stored.label1)
			# if there is only one atom it has to be duplicated for quick_dist2 to work
			if len(atoms1) == 1:
				atoms1=numpy.tile(atoms1, (2,1))
				#print atoms1

			atoms2=numpy.array(stored.label2)
			if len(atoms1) == 1:
				atoms2=numpy.tile(atoms2, (2,1))
				#print atoms2
			
			dist=quick_dist2(self, atoms1, atoms2)
			
			#create pseudoatom at average coordinate of each ensemble and display the distance between them
			avgAtoms1=numpy.average(atoms1,axis=0)
			avgAtoms2=numpy.average(atoms2,axis=0)
			createPseudoatom (avgAtoms1, "tmp_average1", 1)
			createPseudoatom (avgAtoms2, "tmp_average2", 1)
			cmd.distance("avg","tmp_average1 & name ps1","tmp_average2 & name ps1")
			cmd.delete("tmp_average1")
			cmd.delete("tmp_average2")
			
			#cbeta distance if cbeta is present in both selections
			cBeta = []
			if (numberOfAtomsInS1==18 or numberOfAtomsInS1==20) and (numberOfAtomsInS2==18 or numberOfAtomsInS1==20):
				cmd.distance("cBeta", self.residue1_name+" & name CB",self.residue2_name+" & name CB")
				cmd.set("dash_color", "green", "cBeta")
				stored.label1 = []
				stored.label2 = []
				cmd.iterate_state(1, self.residue1_name+" & name CB", 'stored.label1.append((x,y,z))')
				cmd.iterate_state(1, self.residue2_name+" & name CB", 'stored.label2.append((x,y,z))')
				atoms1=numpy.array(stored.label1)
				atoms2=numpy.array(stored.label2)
				cBeta=quick_dist2(self, atoms1, atoms2)
			
			#output of distances and histogram to clipboard or file
			numpy.set_printoptions(threshold=10000000, precision = 2, suppress = True)
			#copy to clipboard
			#create envelope plot
			histogram=numpy.histogram(dist, numpy.arange(100))
			envelopePlot = numpy.zeros((len(dist),2))
			envelopePlot[0:99] = numpy.column_stack((histogram[1][0:len(histogram[1])-1], histogram[0]))
			#put point in mid of bin
			envelopePlot[:,0] += 0.5 
			normEnvelopePlot = numpy.copy(envelopePlot)
			normEnvelopePlot[:,1] = normEnvelopePlot[:,1]/numpy.amax(histogram[0])
			#combine dist and histogram to single array before output
			output=numpy.column_stack((dist, envelopePlot, normEnvelopePlot[:,1]))
			outputStr=numpy.array_str(output)
			outputStr=outputStr.replace("[", "")
			outputStr=outputStr.replace("]", "")
			copyStringToClipboard(outputStr)
			
			if self.writeToFile=='yes':
				numpy.savetxt(self.residue1_name+"-"+self.residue2_name,output, delimiter='\t')
			print calculateStatistics2(dist)
			if len(cBeta) > 0:
				print "Cbeta distance: %3.1f" %cBeta[0]
		
		#Copy and move mode	
        elif self.pick_count == 2 and self.mode == 'Copy & Move':
            print "\n\n\nCopy labels:\n"
            copyAndMove(self, self.residue1_name, self.residue2_name, self.picked_object1, self.numberOfLabel)
                
        #some cleanup after each run
        self.pick_count = 0
        cmd.delete("pk*")
        cmd.delete("sele*")
        cmd.delete("*_tmp*")
        cmd.delete("*tmp*")
        cmd.delete("_indicate*")
        cmd.delete("labelEnvironment*")
        self.cmd.refresh_wizard()
        cmd.set_view(my_view)
        
    
    def delete_last(self):
        try:
            print self.numberOfLabel
            if self.numberOfLabel >= 1:
                cmd.delete(self.object_prefix+str(self.numberOfLabel)+"*")
                self.numberOfLabel-=1
        except pymol.CmdException, pmce:
            print pmce
    
    def toggle_states(self):
        if cmd.get("all_states")=='on':
            self.toggleStatesCaption='Toggle states: OFF'
            cmd.set("all_states",0)
        elif cmd.get("all_states")=='off':
            self.toggleStatesCaption='Toggle states: ON'
            cmd.set("all_states",1)
        self.cmd.refresh_wizard()
            
    def get_panel(self):
        if not ((self.mode == 'Measure') or (self.mode == 'Copy & Move')):
            return [
                    [ 1, 'Mtssl Wizard',''],
                    [ 3, 'Mode: %s'%self.mode,'mode'],
                    [ 3, 'Label: %s'%self.currentLabel,'currentLabel'],
                    [ 3, 'Thoroughness: %s'%self.thoroughness,'thoroughness'],
                    [ 3, 'Cutoff: %3.1f A'%self.cutoff,'cutoff'],
                    [ 3, 'Clashes allowed: %i'%self.clashes,'clashes'],
                    #[ 3, 'Angles to file?: %s'%self.writeToFile,'writeToFile'],
                    [ 2, 'Go!','cmd.get_wizard().run()'],
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
                    
def open_wizard():
    wiz = MtsslWizard()
    cmd.set_wizard(wiz)

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

def copyStringToClipboard(string):
	r = Tk()
	r.withdraw()
	r.clipboard_clear()
	r.clipboard_append(string)
	r.destroy()
	print "Copied to clipboard."

def superpose(residue1,residue2):
    #get the position of the selected residue's O atom
    stored.xyz = []
    cmd.iterate_state(1,residue2+" & name O","stored.xyz.append([x,y,z])")
    #print stored.xyz.pop(0)
    if cmd.pair_fit(residue1+" & name CA",residue2+" & name CA",
                    residue1+" & name N", residue2+" & name N",
                    residue1+" & name C", residue2+" & name C",
                    residue1+" & name CB",residue2+" & name CB"):
        #set the label's O atom to the stored position
        cmd.alter_state(1,residue1+" & name O","(x,y,z)=stored.xyz.pop(0)")
        return True
    else:
        return False

def copyAndMove(self, residue1_name, residue2_name, picked_object1, numberOfLabel):
    label=self.object_prefix+str(numberOfLabel)    
    cmd.copy(label+"_copied", picked_object1)
    superpose(label+"_copied", residue2_name)

def createSelectionMacro(selection):
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
    
def fastMtsslify(self, selected_object, selected_residue, numberOfLabel, thoroughness, cutoff, maxClash):   #generate and check conformations semi-systematically
    #generate the label and superpose onto selected position
    if self.currentLabel == "MTSSL":
    	generateMtssl(numberOfLabel)
    	label="mtssl_"+str(numberOfLabel)
    if self.currentLabel == "PROXYL":
		generateProxyl(numberOfLabel)
		label="proxyl_"+str(numberOfLabel)
		#stored.names = []
		#cmd.iterate_state(1, label, 'stored.names.append(name)')
		#print stored.names
    print label
    print "Attempting superposition..."
    print label
    print selected_residue
    if not superpose(label, selected_residue):
        print "Superposition does not work. Glycine? Mutate to Ala first."
        return
    else:
        print "Superposition worked!"
	print "Trying to find conformations for label %s with thoroughness: %s" % (label, thoroughness)
	#search settings
	numberToFind=100
	if thoroughness == 'painstaking':
		maxNtries=100000
		numberToFind=1000
	elif thoroughness == 'thorough search':
		maxNtries=10000
		numberToFind=200
	elif thoroughness == 'normal search':
		maxNtries=1000
		numberToFind=50
	elif thoroughness == 'quick search':
		maxNtries=50
		numberToFind=50

	#create object with only the atoms around the label to speed everything up 
	protein=selected_object+" &! "+ selected_residue + " within 13.0 of "+label
	#cmd.color('red', protein)
	#print protein
	#cmd.create("labelEnvironment_"+label,protein)
	stored.environmentAtoms = []
	cmd.iterate_state(1, protein, 'stored.environmentAtoms.append((x,y,z))')
	environmentAtoms=numpy.array(stored.environmentAtoms)
	
	#reference atoms are to detect internal clashes
	referenceAtoms=getMovingAtoms(self, label)
	refDist=scipy.spatial.distance.cdist(referenceAtoms, referenceAtoms)
	
	found=0
	ntries=0
	while found < numberToFind and ntries < maxNtries:
		#MTSSL
		#['N', 'O', 'C', 'CA', 'CB', 'SG', 'SD', 'CE', 'C3', 'O1', 'C2', 'N1', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#  0    1    2    3     4     5     6     7     8     9     10    11    12    13    14    15    16    17 
		
		#PROXYL
		#['N', 'O', 'C', 'CA', 'CB', 'SG', 'C1', 'C2', 'N2', 'C3', 'N1', 'O1', 'C10','C4', 'C5', 'C6', 'C7', 'C8', 'C9',  'O3']
		#  0    1    2    3     4     5     6     7     8     9     10    11    12    13    14    15    16    17     18    19 
		
		##Chi1
		#get moving atoms
		movingAtoms=numpy.copy(referenceAtoms)
		axis = numpy.zeros(shape=(2,3))
		#translate moving atoms to origin for easier rotation, will be back-translated at the end
		translationVector=numpy.array(movingAtoms[3])
		backtranslationVector=translationVector
		movingAtoms-=translationVector
		axis[0] = movingAtoms[3]
		axis[1] = movingAtoms[4]
		#rotate moving atoms around axis
		rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
		movingAtoms[5:]=rotatePoints(movingAtoms[5:], rotationMatrix)
		
		##Chi2
		#get axis and translate
		translationVector=numpy.array(movingAtoms[4])
		backtranslationVector+=translationVector
		movingAtoms-=translationVector
		axis[0] = movingAtoms[4]
		axis[1] = movingAtoms[5]
		#rotate moving atoms around axis
		rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
		movingAtoms[6:]=rotatePoints(movingAtoms[6:], rotationMatrix)
		#print movingAtoms
		
		##Chi3
		#get axis
		translationVector=numpy.array(movingAtoms[5])
		backtranslationVector+=translationVector
		movingAtoms-=translationVector
		axis[0] = movingAtoms[5]
		axis[1] = movingAtoms[6]
		#rotate moving atoms around axis
		rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
		movingAtoms[7:]=rotatePoints(movingAtoms[7:], rotationMatrix)
		#print movingAtoms
		
		##Chi4
		#get axis
		translationVector=numpy.array(movingAtoms[6])
		backtranslationVector+=translationVector
		movingAtoms-=translationVector
		axis[0] = movingAtoms[6]
		axis[1] = movingAtoms[7]
		#rotate moving atoms around axis
		rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
		movingAtoms[8:]=rotatePoints(movingAtoms[8:], rotationMatrix)
		#print movingAtoms
		
		##Chi5
		#get axis
		translationVector=numpy.array(movingAtoms[7])
		backtranslationVector+=translationVector
		movingAtoms-=translationVector
		axis[0] = movingAtoms[7]
		axis[1] = movingAtoms[8]
		#rotate moving atoms around axis
		rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
		if self.currentLabel == "MTSSL":
			rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
			movingAtoms[9:]=rotatePoints(movingAtoms[9:], rotationMatrix)
			movingAtoms+=backtranslationVector
		if self.currentLabel == "PROXYL": #this is to account for O3 which is in the rotation axis
			rotationMatrix=setupRotationMatrix(generatePeptideChiAngle(), axis[1])
			#print generatePeptideChiAngle()
			movingAtoms[9:19]=rotatePoints(movingAtoms[9:19], rotationMatrix)
		
		#Chi6 (for PROXYL)
		if self.currentLabel == "PROXYL":
			translationVector=numpy.array(movingAtoms[8])
			backtranslationVector+=translationVector
			movingAtoms-=translationVector
			axis[0] = movingAtoms[8]
			axis[1] = movingAtoms[9]
			#rotate moving atoms around axis
			rotationMatrix=setupRotationMatrix(generateRandomChiAngle(), axis[1])
			movingAtoms[10:19]=rotatePoints(movingAtoms[10:19], rotationMatrix)
			movingAtoms+=backtranslationVector	
		
		if not quickClash(movingAtoms[5:], environmentAtoms, cutoff, maxClash):
			if not internalClash2(movingAtoms, refDist):
				found+=1
				print found,
				#only switch on snuggly fit search for "painstaking"
				if thoroughness == "painstaking":
					vdw = numberOfVdwContacts(self, movingAtoms[5:], environmentAtoms, cutoff)
					thisConformation = [found, vdw]
					self.conformationList.append(thisConformation)	
				createRotamer(self, selected_residue, label, movingAtoms, found)
			else:
				print "i",
		else:
			print ".",
		ntries+=1

	#only switch on snuggly fit search for "painstaking"
	if thoroughness == "painstaking":
		createSnugglyFitConformations(selected_residue, self.conformationList)		
	
	print ""
	print "Found: %i in %i tries." %(found, ntries) 
    print "Done!"
    finalCosmetics(self, label, selected_residue, numberOfLabel,found,thoroughness, maxClash,cutoff)

def calculateStatistics2(distances):
    statisticsResult=""
    #statistics
    average = numpy.average(distances)
    median = numpy.median(distances)
    longest = numpy.amax(distances)
    shortest = numpy.amin(distances)
    
    statisticsResult+= "Average of distribution: %3.2f\n" %average
    statisticsResult+= "Median of distribution: %3.2f\n" %median
    statisticsResult+= "Shortest distance: %3.2f\n" % shortest
    statisticsResult+= "Longest distance: %3.2f" %longest
    return statisticsResult

def numberOfVdwContacts(self, atoms, environmentAtoms, cutoff):
	dist=scipy.spatial.distance.cdist(environmentAtoms, atoms)
	vdwContacts=len(dist[numpy.nonzero((dist > cutoff) & (dist < 4.5 ))])
	return vdwContacts
	

def quick_dist2(self, atoms1, atoms2):
	#take random sample if too many atoms
	if len(atoms1) > 250:
		#atoms1=atoms1[numpy.random.randint(atoms1.shape[0], size = 10),:]
		atoms1=atoms1[numpy.random.permutation(atoms1.shape[0])[:250]]
	if len(atoms2) > 250:
		#atoms2=atoms2[numpy.random.randint(atoms2.shape[0], size = 10),:]
		atoms2=atoms2[numpy.random.permutation(atoms2.shape[0])[:250]]
	dist=scipy.spatial.distance.cdist(atoms1, atoms2)
	dist=numpy.reshape(dist, (-1, 1))
	return dist
	

def createRotamer(self, selected_residue, label, movingAtoms, found):
	if self.currentLabel == "MTSSL":
		names=['N', 'O', 'C', 'CA', 'CB', 'SG', 'SD', 'CE', 'C3', 'O1', 'C2', 'N1', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'] #n1 and c3 swapped, ca and o swapped
		for i in range (0, len(movingAtoms)):
			stored.xyz = []
			stored.xyz = movingAtoms[i]
			cmd.alter_state(1,label +"& name " + names[i],"(x,y,z)=stored.xyz")
		cmd.create(selected_residue, label, 1, found)
		#alter_state 1,selection,(x,y,z)=(newx,newy,newz)
	if self.currentLabel == "PROXYL":
		names=['N', 'O', 'C', 'CA', 'CB', 'SG', 'C1', 'C2', 'N2', 'C3', 'N1', 'O1', 'C10','C4', 'C5', 'C6', 'C7', 'C8', 'C9',  'O3'] #swapped
		for i in range (0, len(movingAtoms)):
			stored.xyz = []
			stored.xyz = movingAtoms[i]
			cmd.alter_state(1,label +"& name " + names[i],"(x,y,z)=stored.xyz")
		cmd.create(selected_residue, label, 1, found)
	
def internalClash2(atoms, refDist):
	#distances in new rotamer
	dist=scipy.spatial.distance.cdist(atoms, atoms)
	#create Boolean array with elements that describe if a distance changes or not
	changingDistances = numpy.absolute(numpy.round(numpy.subtract(dist,refDist),2)) > 0
	#multiply by Boolean area to make all constant distances zero
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

def getMovingAtoms(self, selection):
	if self.currentLabel == "MTSSL":
		#This is the order of atoms when iterate is used on mtssl in pymol:
		#['N', 'CA', 'C', 'O', 'CB', 'SG', 'SD', 'CE', 'N1', 'O1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#  0    1    2    3     4     5     6     7     8     9     10    11    12    13    14    15    16    17  
		#This is the order we want
		#['N', 'O', 'C', 'CA', 'CB', 'SG', 'SD', 'CE', 'C3', 'O1', 'C2', 'N1', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
		#  0    1    2    3     4     5     6     7     8     9     10    11    12    13    14    15    16    17  
		cmd.select("tmp", selection)
		stored.movingAtoms = []
		cmd.iterate_state(1, 'tmp', 'stored.movingAtoms.append((x,y,z))')
		#swap N1 with C3 and CA with O, so that all atoms which serve as rotation axes are in sequence
		n1=stored.movingAtoms[8]
		c3=stored.movingAtoms[11]
		
		ca=stored.movingAtoms[1]
		o=stored.movingAtoms[3]
		
		stored.movingAtoms[8]=c3
		stored.movingAtoms[11]=n1
		
		stored.movingAtoms[1]=o
		stored.movingAtoms[3]=ca
		
		movingAtoms=numpy.array(stored.movingAtoms)

	if self.currentLabel == "PROXYL":
		#This is the order of atoms when iterate is used on proxyl in pymol:
		#['N', 'CA', 'C', 'O', 'CB', 'SG', 'C1', 'N1', 'O1', 'C2', 'N2', 'C3', 'O3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10']
		#  0    1     2    3    4     5     6     7     8     9     10    11    12    13    14    15    16    17     18    19 
		#This is the order we want:
		#['N', 'O', 'C', 'CA', 'CB', 'SG', 'C1', 'C2', 'N2', 'C3', 'N1', 'O1', 'C10','C4', 'C5', 'C6', 'C7', 'C8', 'C9',  'O3']
		#  0    1    2    3     4     5     6     7     8     9     10    11    12    13    14    15    16    17     18    19 
		
		cmd.select("tmp", selection)
		stored.movingAtoms = []
		cmd.iterate_state(1, 'tmp', 'stored.movingAtoms.append((x,y,z))')
		#swap N1 with C3 and CA with O, so that all atoms which serve as rotation axes are in sequence
		
		ca=stored.movingAtoms[1]
		o=stored.movingAtoms[3]
		n1=stored.movingAtoms[7]
		o1=stored.movingAtoms[8]
		c2=stored.movingAtoms[9]
		n2=stored.movingAtoms[10]
		c3=stored.movingAtoms[11]
		o3=stored.movingAtoms[12]
		c10=stored.movingAtoms[19]
		
		stored.movingAtoms[1]=o
		stored.movingAtoms[3]=ca
		stored.movingAtoms[7]=c2
		stored.movingAtoms[8]=n2
		stored.movingAtoms[9]=c3
		stored.movingAtoms[10]=n1
		stored.movingAtoms[11]=o1
		stored.movingAtoms[12]=c10
		stored.movingAtoms[19]=o3
		
		movingAtoms=numpy.array(stored.movingAtoms)
	return movingAtoms

def createPseudoatom (coordinates, objectName, state):
	x=float(coordinates[0])
	y=float(coordinates[1])
	z=float(coordinates[2])
	posString="[%3.2f,%3.2f,%3.2f]" % (x,y,z)
	cmd.pseudoatom(pos=posString, object=objectName, state=state)

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

def createSnugglyFitConformations(selected_residue, conformationList): #select conformations with the best fit to the molecular surface, rank them and create an object
    print "Snuggliest fit(s):"
    #calculate average atom count of all conformations
    atomCountSum=0
    snugglyFitList=[]
    bestSnugglyFitAtomCount=0
    for x in range (0,len(conformationList)):
        thisConformation=conformationList[x]
        atomCountSum+=thisConformation[1]
        if thisConformation[1] > bestSnugglyFitAtomCount:
            bestSnugglyFitAtomCount=thisConformation[1]
    averageAtomCount=atomCountSum/len(conformationList)
    
    #generate snugglyFitList: take only those conformations whose atom count is > 0.75 of top peak and higher than the average
    counter=1
    for x in range (0,len(conformationList)):
        thisConformation=conformationList[x]
        if thisConformation[1] > 0.75 * bestSnugglyFitAtomCount and thisConformation[1] > averageAtomCount:
            snugglyFitList.append({'vdw':thisConformation[1],
                                   'state':thisConformation[0]})
            cmd.create(selected_residue+"_snuggly", selected_residue, thisConformation[0], counter)
            counter+=1
            
    #sort snugglyFitList so that best fitting conformation is on top   
    snugglyFitList = sorted(snugglyFitList, key=itemgetter('vdw'))
    snugglyFitList.reverse()
    #print out the list
    if len(snugglyFitList)>0:
        for i in range (0,len(snugglyFitList)):
            print "Conformation %i: \t\t%i \t\t\t vdW contacts" %(snugglyFitList[i]['state'],snugglyFitList[i]['vdw'])
    print "Average vdW contacts of all possible conformations: ",averageAtomCount

def finalCosmetics(self, label, selected_residue,numberOfLabel, found, thoroughness,maxClash,cutoff): #make everything look nice
	if found >= 1:
		print "Found %i conformations." %found
		#show all conformations and do some coloring
		cmd.set("all_states",1)
		self.toggleStatesCaption='Toggle states: ON'
		cmd.color("blue",selected_residue)
		cmd.color("red", "name O1")
		cmd.disable(label)
		cmd.disable(label+"_internal_clash")
		cmd.show("spheres", "name O1")
		cmd.set("sphere_scale", "0.2")
		identifierLabel="%s|%s|%1.2f|%i|%s" %(selected_residue, self.currentLabel, cutoff, maxClash, thoroughness)
		print identifierLabel
		#mark average N1 position with pseudoatom
		stored.label = []
		cmd.iterate_state(0, selected_residue+" & name N1", 'stored.label.append((x,y,z))')
		atoms1=numpy.array(stored.label)
		#create pseudoatom at average coordinate of each ensemble
		avgAtoms=numpy.average(atoms1,axis=0)
		createPseudoatom (avgAtoms, selected_residue+"_label", 1)
		cmd.set("sphere_scale", "0.5", selected_residue+"_label")
		cmd.label(selected_residue+"_label", `identifierLabel`)
		cmd.show("label")
		cmd.show("spheres", "name PS1")
		print label+"*"+","+"labelEnvironment_"+label
		cmd.delete(label+"*")
		cmd.group(self.object_prefix+str(numberOfLabel), label+"*"+","+"labelEnvironment_"+label+","+selected_residue+"*")
	else:
		print "Sorry, I did not find anything. Your options are:\n1) Try again or\n2) Try to increase 'thoroughness' (now: "+str(thoroughness)+"), \nincrease'maxClash' (now: "+str(maxClash)+") or \ndecrease cutoff (now: "+str(cutoff)+")."
		cmd.delete(label+"*")
		print label

def generateProxyl(numberOfLabel): #create a mtssl label object.
    cmd.read_pdbstr("""HEADER    PROXYL\n
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
ATOM     38 C10  IA1 A   1      -7.356  -3.824   5.436  1.00 20.00           C\n""", "proxyl_" + str(numberOfLabel))  

def generateMtssl(numberOfLabel): #create a mtssl label object.
    cmd.read_pdbstr("""HEADER    MTSSL\n
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
ATOM     36 O    R1A A   1      -0.298   2.670  -0.967  1.00 20.00           O\n""", "mtssl_" + str(numberOfLabel))