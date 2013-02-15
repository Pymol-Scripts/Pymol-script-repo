"""
----------------------------------------------------------------------
----------------------------------------------------------------------
--- mtsslDock: dock proteins based on PELDOR distance constraints  --- 
Author	: Gregor Hagelueken
Date	: Feb 2013
Version : 1.0
Mail	: hagelueken'at'pc.uni-bonn.de
----------------------------------------------------------------------
----------------------------------------------------------------------
"""

import pymol
import numpy
import scipy.spatial.distance
import random, time, math
import os
from pymol import cmd
from pymol import util
from pymol import stored
from operator import itemgetter

def createPseudoatom (coordinates, objectName, state):
		x=float(coordinates[0])
		y=float(coordinates[1])
		z=float(coordinates[2])
		posString="[%3.2f,%3.2f,%3.2f]" % (x,y,z)
		cmd.pseudoatom(pos=posString, object=objectName, state=state)

def generateChromosome():
	#translation between -50 and 50
	translation = 50 * numpy.random.uniform(-1, 1, size=3)
	rotation=360 * numpy.random.uniform(0, 1, size=3)
	fitness=numpy.zeros(1)
	#print numpy.concatenate((rotation, translation), axis = 1)
	return numpy.concatenate((rotation, translation, fitness), axis = 1)

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
	#print angle
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

def getDistances(proteinA, proteinB):
	distances = [numpy.linalg.norm(proteinA[0]-proteinB[0]),
			 	 numpy.linalg.norm(proteinA[0]-proteinB[1]),
			 	 numpy.linalg.norm(proteinA[0]-proteinB[2]),
			 	 numpy.linalg.norm(proteinA[1]-proteinB[0]),
			 	 numpy.linalg.norm(proteinA[1]-proteinB[1]),
			 	 numpy.linalg.norm(proteinA[1]-proteinB[2]),
			 	 numpy.linalg.norm(proteinA[2]-proteinB[0]),
			 	 numpy.linalg.norm(proteinA[2]-proteinB[1]),
			 	 numpy.linalg.norm(proteinA[2]-proteinB[2])]
	return distances

def quickClash(atomsProteinA, atomsProteinB, cutoff):
	#atomsProteinA=atomsProteinA[numpy.random.permutation(atomsProteinA.shape[0])[:11]]
	#atomsProteinB=atomsProteinB[numpy.random.permutation(atomsProteinB.shape[0])[:11]]
	#print len(atomsProteinA), len(atomsProteinB)
	dist=scipy.spatial.distance.cdist(atomsProteinA, atomsProteinB)
	#print dist
	clashes=len(numpy.nonzero(dist < cutoff)[0])
	#print clashes
	return clashes

def getFitness(expDistances, trialDistances, atomsProteinA, atomsProteinB, checkDistances):
	diff = 0
	clashes = 0
	for i in range (0, len(expDistances)):
		diff += numpy.absolute(expDistances[i]-trialDistances[i])
	if checkDistances:
		clashes = quickClash(atomsProteinA, atomsProteinB, 2.0)	
	diff = diff + clashes
	#print diff, clashes
	return diff

def getGeneticOperator():
	if random.choice([True, False]):
		if random.choice([True, False]):
			return "smallCreepMutation"
		else:
			return "randomMutation"
	else:
		if random.choice([True, False]):
			return "singlePointCrossover"
		else:
			return "exchangeCrossover"

def smallCreepMutation(parent):
	childs = []
	position = random.randrange(6)
	oldValue = parent[position]
	newValue = 0
	#rotation
	if position < 3:
		creep = 0
		while creep == 0:
			creep = random.randrange(-5.0, 5.0)
		newValue = oldValue + creep
		if newValue > 360:
			newValue -= 360
		elif newValue < 0:
			newValue += 360
	#translation
	else:
		creep = 0
		while creep == 0:
			creep = random.randrange(-2.0, 2.0)
		newValue = oldValue + creep
	parent[position] = newValue
	childs.append(parent)
	return childs

def randomMutation(parent):
	childs = []
	position = random.randrange(6)
	oldValue = parent[position]
	newValue = 0
	#rotation
	if position < 3:
		newValue = random.randrange(360)
	#translation
	else:
		newValue = random.randrange(-50, 50)
	parent[position] = newValue
	childs.append(parent)
	return childs
	
def singlePointCrossover(parent1, parent2):
	childs = []
	position = random.randrange(6)
	child1 = generateChromosome()
	child2 = generateChromosome()
	for i in range (0, 5):
		if i <= position:
			child1[i]=parent1[i]
			child2[i]=parent2[i]
		else:
			child1[i]=parent2[i]
			child2[i]=parent1[i]
	childs.append(child1)
	childs.append(child2)
	return childs
	
def exchangeCrossover(parent1, parent2):
	childs = []
	position = random.randrange(6)
	child1 = numpy.copy(parent1)
	child2 = numpy.copy(parent2)
	tmpValue = child1[position]
	child1[position] = child2[position]
	child2[position] = tmpValue
	child1[6] = 0
	child2[6] = 0
	return childs

def createTrialPosition(proteinB, chromosome):
	#rotate around x
	trialProteinB = numpy.copy(proteinB)
	#print chromosome[0]
	rotationMatrix = setupRotationMatrix(chromosome[0], numpy.array([1,0,0]))
	trialProteinB = rotatePoints(trialProteinB, rotationMatrix)
	#rotate around y
	rotationMatrix = setupRotationMatrix(chromosome[1], numpy.array([0,1,0]))
	trialProteinB = rotatePoints(trialProteinB, rotationMatrix)
	#rotate around z
	rotationMatrix = setupRotationMatrix(chromosome[2], numpy.array([0,0,1]))
	trialProteinB = rotatePoints(trialProteinB, rotationMatrix)
	
	#transform along x, y, z
	#print trialProteinB
	trialProteinB = trialProteinB + numpy.array([chromosome[3],0,0])
	trialProteinB = trialProteinB + numpy.array([0,chromosome[4],0])
	trialProteinB = trialProteinB + numpy.array([0,0,chromosome[5]])
	return trialProteinB

def moveProtein(proteinBPymolString, solutionPymolString, solution, state):
	#recreate protein B at solution coordinates
	cmd.create(solutionPymolString, proteinBPymolString, 1, state)
	#translate to origin with proteinBcog. IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system!
	cmd.translate(list(-1*proteinBcog.reshape(-1,)), solutionPymolString,state,0,None)
	
	#rotate and translate according to solution
	translation = solution[3:6]
	#print list(translation.reshape(-1,))
	rotX = solution[0]
	#IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system! Also set rotation origin to 0,0,0!
	cmd.rotate([1, 0, 0], rotX, solutionPymolString,state,0,None,[0,0,0])
	rotY = solution[1]
	cmd.rotate([0, 1, 0], rotY, solutionPymolString,state,0,None,[0,0,0])
	rotZ = solution[2]
	cmd.rotate([0, 0, 1], rotZ, solutionPymolString,state,0,None,[0,0,0])
	cmd.translate(list(translation.reshape(-1,)), solutionPymolString,state,0,None)	

def printChromosomes(chromosomes):
	if len(chromosomes) == 1:
		print "rotX: %1.2f\t rotY: %1.2f\t rotZ: %1.2f\t tX: %1.2f\t tY: %1.2f\t tZ: %1.2f\t Fit:%1.2f" %(chromosomes[0],chromosomes[1],chromosomes[2],chromosomes[3],chromosomes[4],chromosomes[5], chromosomes[6])
	else:
		for chromosome in chromosomes:
			print "rotX: %1.2f\t rotY: %1.2f\t rotZ: %1.2f\t tX: %1.2f\t tY: %1.2f\t tZ: %1.2f\t Fit:%1.2f" %(chromosome[0],chromosome[1],chromosome[2],chromosome[3],chromosome[4],chromosome[5], chromosome[6])

def evolve(generations, chromosomes, mutationFrequency,rigidBody, scoreClashes):
	for i in range (0,generations):
		print "%i of %i" %(i, generations)
		#rank
		sortedChromosomes = sorted(chromosomes, key=itemgetter(6))
		#remove 5% unfittest and produce new offspring
		replaceFromIndex = int(len(sortedChromosomes)*mutationFrequency)
		survivingChromosomes = sortedChromosomes[0:len(sortedChromosomes)-replaceFromIndex]
		for j in range (0, replaceFromIndex):
			geneticOperator = ""
			if not rigidBody:
				geneticOperator = getGeneticOperator()
			else:
				while not geneticOperator == "smallCreepMutation":
					geneticOperator = getGeneticOperator()
					#print survivingChromosomes
			childs = []
			if geneticOperator == "smallCreepMutation":
				if rigidBody:
					#parent = survivingChromosomes[0]
					parent = survivingChromosomes[random.randrange(3)]
					print parent
				else:
					parent = survivingChromosomes[random.randrange(len(survivingChromosomes))]
				childs = smallCreepMutation(numpy.copy(parent))
			elif geneticOperator == "randomMutation":
				parent = survivingChromosomes[random.randrange(len(survivingChromosomes))]
				childs = randomMutation(numpy.copy(parent))
			elif geneticOperator == "singlePointCrossover":
				parent1 = survivingChromosomes[random.randrange(len(survivingChromosomes))]
				parent2 = survivingChromosomes[random.randrange(len(survivingChromosomes))]
				childs = singlePointCrossover(numpy.copy(parent1), numpy.copy(parent2))
			elif geneticOperator == "exchangeCrossover":
				parent1 = survivingChromosomes[random.randrange(len(survivingChromosomes))]
				parent2 = survivingChromosomes[random.randrange(len(survivingChromosomes))]
				childs = exchangeCrossover(numpy.copy(parent1), numpy.copy(parent2))
			#Determine fitness of childs
			for child in childs:
				trialProteinB = createTrialPosition(proteinB, child)
				trialDistances = getDistances(proteinA, trialProteinB)
				if scoreClashes:
					trialTmpBatoms = createTrialPosition(numpy.copy(tmpBatoms), child)
					fitness = getFitness(expDistances, trialDistances, tmpAatoms, trialTmpBatoms, True)
				else:
					fitness = getFitness(expDistances, trialDistances, tmpAatoms, tmpBatoms, False)
				child[6] = fitness
			survivingChromosomes+=childs
			survivingChromosomes = sorted(survivingChromosomes, key=itemgetter(6))
			if rigidBody:
				printChromosomes(survivingChromosomes)
		#prepare for next generation
		chromosomes = survivingChromosomes
	return chromosomes
		

numberOfChromosomes = 400
myView = cmd.get_view()
#RdxR label cogs
proteinA = numpy.array([[-37.88999938964844, -32.619998931884766, -1.9600000381469727], [-14.800000190734863, -1.6299999952316284, 10.5600004196167], [21.860000610351562, -26.639999389648438, 16.459999084472656]])
proteinAcog = numpy.average(proteinA,axis=0)

#Rdx label cogs
proteinB = numpy.array([[-12.460000038146973, -20.09000015258789, 24.219999313354492], [-26.829999923706055, -42.65999984741211, 19.31999969482422], [0.10999999940395355, -35.5099983215332, 20.15999984741211]])
proteinBcog = numpy.average(proteinB,axis=0)

expDistances = getDistances(proteinA, proteinB)
#print expDistances

#move both to origin
proteinA = proteinA - proteinAcog
#for point in proteinA:
#	createPseudoatom(point, "ProteinAinitial", 1)

proteinB = proteinB - proteinBcog
#for point in proteinB:
#	createPseudoatom(point, "ProteinBinitial", 1)


cmd.create("tmpA", "2v3b & chain A", 1, 1)
cmd.translate(list(-1*proteinAcog.reshape(-1,)), "tmpA",1,0,None)

stored.tmpAatoms = []
cmd.iterate_state(1, "tmpA", 'stored.tmpAatoms.append((x,y,z))')
tmpAatoms = numpy.array(stored.tmpAatoms)

cmd.create("tmpB", "2v3b & chain B", 1, 1)
cmd.translate(list(-1*proteinBcog.reshape(-1,)), "tmpB",1,0,None)

stored.tmpBatoms = []
cmd.iterate_state(1, "tmpB", 'stored.tmpBatoms.append((x,y,z))')
tmpBatoms = numpy.array(stored.tmpBatoms)

chromosomes = []

for i in range (0,numberOfChromosomes):
	chromosomes.append(generateChromosome())

#calculate fitness of initial population
i=0
for chromosome in chromosomes:
	print "Initial scoring..."
	trialProteinB = createTrialPosition(proteinB, chromosome)
	trialDistances = getDistances(proteinA, trialProteinB)
	fitness = getFitness(expDistances, trialDistances, tmpAatoms, tmpBatoms, False)
	#print fitness
	chromosome[6] = fitness
	i+=1

generation1 = evolve(1000, chromosomes, 0.05, False, False)
generation1 = sorted(generation1, key=itemgetter(6))

#select 10% fittest, score conformations including clashes and respawn population then do rigid body
generation1 = generation1[0:int(len(generation1)*0.1)]
print generation1
i = 0
for chromosome in generation1:
	print "Adding clashes to fitness function..."
	print "%i of %i" %(i, len(generation1))
	trialProteinB = createTrialPosition(proteinB, chromosome)
	trialDistances = getDistances(proteinA, trialProteinB)
	trialTmpBatoms = createTrialPosition(numpy.copy(tmpBatoms), chromosome)
	fitness = getFitness(expDistances, trialDistances, tmpAatoms, trialTmpBatoms, True)
	#print fitness
	chromosome[6] = fitness
	i+=1
print "Final rigid body refinement..."
generation2 = evolve(20, generation1, 0.2, True, True)
generation2 = sorted(generation2, key=itemgetter(6))
solutions = generation2[0:10]
counter = 1
print ""
printChromosomes(solutions)
for solution in solutions:
	solutionProteinB = createTrialPosition(proteinB, solution)
	#recreate protein B at solution coordinates
	moveProtein("2v3b & chain B", "solution", solution, counter)
	#translate back to match with protein A original position
	cmd.translate(list(proteinAcog.reshape(-1,)), "solution",counter,0,None)
	counter+=1

#also recreate protein A at origin.
#cmd.create("solution_%i_proteinA" %counter, "2v3b & chain A", 1, 1)
#cmd.translate(list(-1*proteinAcog.reshape(-1,)), "solution_%i_proteinA" %counter,1,0,None)
cmd.set("all_states", 1)
cmd.delete("tmp*")
cmd.set_view(myView)