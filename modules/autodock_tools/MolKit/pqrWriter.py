#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

from MolKit.moleculeWriter import MoleculeWriter
from MolKit.protein import Protein, Chain, ChainSet
from MolKit.molecule import Atom
from MolKit.tree import TreeNode, TreeNodeSet
import os

class PqrWriter(MoleculeWriter):
    """Class to write data records from a molecule tree to a prq file."""
    

    def write(self, fileName, node, sort=None):
	    """Takes a filename and TreeNode or TreeNodeSet instance.  If no
	    filename extension is provided, a '.pqr' extension is added. For
	    the node or set, the records for the whole protein are written,
	    but ATOM records are written only for those atoms contained
	    within and below that node or set.
	    if a sort function is specified, the list of nodes is sorted
	    using this function. This function has to return -1, 0, 1 if
	    the first argument is smaller, equal or larger then the second
	    argument"""

	    self.numCoord = 0
	    filename = fileName
	    root_ext = os.path.splitext('%s' %filename)
	    if root_ext[1]=='':
		    filename = '%s.pqr' %fileName
	    assert isinstance(node, TreeNode) or isinstance(node, TreeNodeSet)
       	    if isinstance(node, TreeNode):
	        mol = node.top   
	    elif isinstance(node, TreeNodeSet):
	        mol = node.top.uniq()[0]
	    else:
	        return

	    file = open('%s' % filename, 'w')
	    if sort:
	        node.sort(sort)
	    self.write_AtomSet(file, node)
	    file.close()
    


    def write_atom(self, f, atm):
        """Takes a file object and an Atom instance.  Writes the atom
	record to the file."""
	    
	spaceStr = "          "
	#columns1-6
	if atm.hetatm==0: f.write('ATOM  ')
	else: f.write('HETATM')
	#columns 7-11 + A SPACE
	f.write('%5i ' % atm.number)
	#columns 13-16 
	spaceChar = None

	if len(atm.name)==4:
	    if atm.element=='H':
	        #name = atm.name[-1]+atm.name[1:]
	        name = atm.name[-1]+atm.name[:-1]
	        f.write('%-4.4s' % name)
            else:
	        name = atm.name
	        f.write('%4.4s' % name)
	elif len(atm.element)==2: #else of  if(atm.name)==4
	        f.write('%-4.4s' % atm.element)
	else:
	    f.write(' %-3s' % atm.name)

	#columns 17
	if not spaceChar:
	    f.write('%1.1s' %spaceStr)
	else:
	    f.write('%1.1s' %spaceChar[-1])

	#columns 18-20 SPACE 22
	if hasattr(atm, 'parent') and hasattr(atm.parent, 'type'):
	    f.write('%3.3s ' %atm.parent.type)
            # no chain ID
            f.write('%1.1s' %spaceStr)
## 	    if hasattr(atm.parent, 'parent') and hasattr(atm.parent.parent, 'id'):
## 	        f.write('%1.1s' %atm.parent.parent.id)
## 	    else:
## 	        f.write('%1.1s' %spaceStr)
	else:
	    f.write('%5.5s' %spaceStr)

	#columns 23-26
	if hasattr(atm.parent, 'number'):
	    f.write('%4.4s' % atm.parent.number)
	else:
	    f.write('%4.4s' % spaceStr)
	#columns 27 plus 3 SPACES
	if hasattr(atm.parent, 'icode'):
	    z=atm.parent.icode+'   '
        else: 
	    z="    "
	f.write('%4.4s' %z)
	#columns 31-38,39-46, 47-54
	for i in atm.coords:
	    f.write('%8.3f' %i)

	f.write('%8.3f' %atm.charge)
        if hasattr(atm, 'pqrRadius'):
            f.write('%8.3f' %atm.pqrRadius)
        else:
            f.write('%8.3f' %atm.radius)

	f.write('\n')


    def write_AtomSet(self, file, node):
        """Takes a file object and a TreeNode or TreeNodeSet instance.
	For each Atom in node, write_AtomSet calls method write_atom.
	If node is a Chain or higher, write_AtomSet calls write_TER in
	between chains or at the end of a chain."""

	if isinstance(node, Protein) or isinstance(node,Chain) or \
	   isinstance(node,ChainSet):

	    chains = node.findType(Chain)
	    for c in chains:
	        atmset = c.findType(Atom)
		try:
		    atmset.charge
		except:
		    print 'ERROR: atoms with missing charge found'
		    return 'Error'

		try:
		    atmset.radius
		except:
		    print 'ERROR: atoms with missing radius found'
		    return 'Error'

		for a in atmset:
		    self.write_atom(file, a)
		    self.numCoord = self.numCoord + 1
	        self.write_TER(file, atmset[-1])
		self.numTer = self.numTer + 1
	else:
	    atmset = node.findType(Atom)
	    try:
	        atmset.charge
	    except:
	        print 'ERROR: atoms with missing charge found'
		return 'Error'
	    try:
	        atmset.radius
	    except:
	        print 'ERROR: atoms with missing radius found'
		return 'Error'

	    for a in atmset:
	        self.write_atom(file, a)
		self.numCoord = self.numCoord + 1	    
