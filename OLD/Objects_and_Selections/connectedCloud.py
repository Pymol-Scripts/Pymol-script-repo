# -*- coding: utf-8 -*-
def connectedCloud( origSel, subSel, radius=1, doPrint=True ):
	"""
	connectedCloud -- find all atoms, by object in subSel connected to
		origSel by radius Ang.  In other words, start from origSel
		and find all objects <=radius Ang. from it & add that to the
		result.  Once nothing more is added, you have your cloud of
		points.
 
	origSel,
		the selection around which we start expanding.
 
	subSel,
		the subset of objects to consider
 
	radius,
		the maximum radius in which to consider; anything >radius
		Ang. from any objet in the current selection is ignored (unless
		it is added later by some other path).
 
	doPrint
		print some values before returning the selection name
 
	Author:
		Jason Vertrees, 2009.
 
	"""
	cName = "__tempCloud__"
	cmd.select(cName, origSel)
 
	# we keep track of cloud size based on the number of atoms in the selection; oldCount
	# is the previous iteration and curCount is after we add everything for this iteration.
	oldCount=0
	curCount=cmd.count_atoms( cName )
 
	# flag to stop at 1000 iterations. if something goes wrong we don't want to hijack the computer.
	flag=0
 
	# while something was added in this iteration
	while oldCount != curCount:
		flag += 1
		# grow current cloud
		cmd.select( cName, "bo. " + subSel + " within 1 of " + cName )
		# update atom counts in the cloud
		oldCount=curCount
		curCount=cmd.count_atoms(cName)
 
		if ( flag > 1000 ):
			print "Warning: I gave up.  I grew the cloud 1000 times and it kept growing."
			print "Warning: If this is valid (weird), then increases the 'flag' in the source"
			print "Warning: code to something >1000."
			break
 
	# create the object for getting its surface area
	cmd.create("connectedCloud", cName)
 
	# print cloud info
	if doPrint:
		print "Number of atoms in cloud:", str(cmd.count_atoms(cName))
		print "Surface area of cloud:", str(cmd.get_area("connectedCloud"))
 
	# delete our work 
	cmd.delete("connectedCloud")
	# rename and return
	cmd.set_name(cName, "connectedCloud")
	return "connectedCloud"
 
cmd.extend("connectedCloud", connectedCloud)

