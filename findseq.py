"""
DESCRIPTION:
Given a sequence/regex to find, select those
matching amino acids in the protein.

USAGE:
findseq needle, haystack[, selName[, het[, firstOnly]]]

PARAMS:
needle (string)
		the sequence of amino acids to match and select
		in the haystack.  This can be a sequence of amino
		acids, or a string-style regular expression.  See
		examples.

hastack (string or PyMOL selection)
		name of the PyMOL object/selection in which 
		to find the needle.

selName (string; defaults to None)
		This is the name of the selection to return.  If selName
		is left blank (None), then the selection name will be
		foundSeqXYZ where XYZ is some random number; if selName is
		"sele" the usual PyMOL "(sele)" will be used; and, lastly,
		if selName is anything else, that name will be used verbatim.

het (0 or 1; defaults to 0)
		This boolean flag allows (1) or disallows (0) heteroatoms
		from being considered.

firstOnly (0 or 1; defaults to 0)
		Subsequences or motifs might be repeated, this controls how we
		consider multiple matches.  If firstOnly is False (0) then we return
		all found subsequences; if firstOnly is True (1), then we just return
		the first found sequence.

RETURNS:
a newly created selection with the atoms you sought.  If there are
more than two contiguous regions, then a newly created group is
returned with each contiguous segment its own selection.

EXAMPLE:
# find SPVI in 1h12, foundSeqXYZ as return name
findseq SPVI, 1h12

# find FATEW and make it (sele).
findseq FATEW, 1g01, sele

# find the regular expression GMS.*QWY in 1a3h
# and put the return value in (sele).
fetch 1a3h
# this ends up finding the sequence, GMSSHGLQWY
findseq GMS.*QWY, 1a3h, sele

NOTES:
Assumes we're using the ONE LETTER amino acid abbreviations.

AUTHOR:
Jason Vertrees, 2009.
"""

from pymol import cmd
import re,types,random

def findseq(needle, haystack, selName=None, het=0, firstOnly=0):
	# set the name of the selection to return.
	if selName == None:
		rSelName = "foundSeq" + str(random.randint(0,32000))
		selName = rSelName
	elif selName == "sele":
		rSelName = "sele"
	else:
		rSelName = selName

	# input checking
	if not checkParams(needle, haystack, selName, het, firstOnly):
		print "There was an error with a parameter.  Please see"
		print "the above error message for how to fix it."
		return None

	one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
	'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
	'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
	'GLY':'G', 'PRO':'P', 'CYS':'C'}

	# remove hetero atoms (waters/ligands/etc) from consideration?
	if het:
		cmd.select( "__h", "br. " + haystack )
	else:
		cmd.select( "__h", "br. " + haystack + " and not het" )

	# get the AAs in the haystack
	aaDict = { 'aaList' : [] }
	cmd.iterate("(name ca) and __h","aaList.append((resi,resn))",space=aaDict)

	IDs = map( lambda x: int(x[0]), aaDict['aaList'] )
	AAs = ''.join(map( lambda x: one_letter[x[1]], aaDict['aaList'] ))

	reNeedle = re.compile( needle.upper() )
	it = reNeedle.finditer( AAs )

	# make an empty selection to which we add residues
	cmd.select( rSelName, 'None')

	for i in it:
		(start,stop) = i.span()
		cmd.select( rSelName, rSelName + " __h and i. " + str(IDs[start]) + "-" + str(IDs[stop-1]))
		if int(firstOnly):
			break
	cmd.delete("__h")
	return rSelName
cmd.extend("findseq", findseq )

def checkParams(needle,haystack,selName,het,firstOnly):
	"""
	This is just a helper function for checking the user input
	"""
	# check Needle
	if len(needle)==0 or type(needle)!=types.StringType:
		print "Error: Please provide a string 'needle' to search for."
		print "Error: For help type 'help motifFinder'."
		return False

	# check Haystack
	if len(haystack)==0 or type(haystack)!=types.StringType:
		print "Error: Please provide valid PyMOL object or selection name"
		print "Error: in which to search."
		print "Error: For help type 'help motifFinder'."
		return False

	# check het
	try:
		het = bool(int(het))
	except ValueError:
		print "Error: The 'het' parameter was not 0 or 1."
		return False

	# check first Only
	try:
		firstOnly = bool(int(het))
	except ValueError:
		print "Error: The 'firstOnly' parameter was not 0 or 1."
		return False

	# check selName
	if type(selName)!=types.StringType:
		print "Error: selName was not a string."
		return False
	return True
