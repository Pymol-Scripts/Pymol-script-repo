#
#
#
#
#
#
#$Id: groParser.py,v 1.11 2006/12/15 19:53:29 annao Exp $
# 

"""
Module Gromacs Parser.
"""


import os,types    
from MolKit.moleculeParser import MoleculeParser
from MolKit.protein import Protein, Chain, ChainSet, Residue, ResidueSet, ProteinSet
from MolKit.molecule import Atom, AtomSet, Bond, BondSet, HydrogenBond
from PyBabel.babelElements import babel_elements


class groParser( MoleculeParser ):

    def __init__( self, filename=None, allLines=None ):
        """Constructor for groParser: adopted form PdbParser"""
        MoleculeParser.__init__( self, filename, allLines )       
        
    def parse( self, objClass=Protein ):
        if self.allLines is None and self.filename:
            self.readFile()
            if self.allLines is None or len(self.allLines)==0:
                return

        mol = Protein()
        self.mol = mol
        molList = mol.setClass()
        molList.append( mol )
        current_residue_number = None
        current_chain = None
        current_residue = None
        number_of_atoms = int(self.allLines[1][:5])

        self.configureProgressBar( init=1, mode='increment', 
                                  authtext='parse atoms', max=number_of_atoms )
        
        
        current_chain = Chain( id='GRO',)
        #FIX this: The existence of allAtoms attribute (and the fact that it is an empty set rather than all atoms in the chain) causes getNodesByMolecule() to return wrong values
        if hasattr(current_chain, "allAtoms"):
            del(current_chain.allAtoms)
        #current_chain = Chain( id='GRO',parent = mol)
        mol.adopt( current_chain, setChildrenTop=1 )
         
        for index in range( 2,number_of_atoms+2 ):              
            residue_number = int(self.allLines[index][:5])
            if residue_number!=current_residue_number:# 
                #current_chain should adopt the current residue if there is one
                #create new residue
                res_type = self.allLines[index][5:10]
                residue_type = res_type.split(' ')[0]
                
                current_residue = Residue( type=residue_type, number=residue_number )
                current_residue_number = residue_number
                if current_residue is not None:    #REMEMBER TO ADOPT THE LAST ONE!!!
                     
                    current_chain.adopt( current_residue, setChildrenTop=1 )
                            
            n = self.allLines[index][10:15]
            name = n.split(' ')[-1]
            element = name 
            
            if element in babel_elements.keys():
                element = element

            else:
                 
                if residue_type == "System" or residue_type == "SOL":                 
                    #if element[1] == 'W':
                    #          element = 'H'
                        #   group is treated as one particle
                    #else:
                    element = element[0]

                elif element[:2] == 'Me':
                    element = 'C'
                else:
                    element = element[0]
                
            #if len(element)>1:
            #    if type(element[1]) == types.StringType:
            #        
            #        if element[1] == element[1].lower():
            #            element =element
            #        else:
            #            element = element[0]
            #            
            #    else:    
            #        element = element[0]
                
            atom = Atom( name, current_residue, element, top=mol )
            c =  self.allLines[index][15:20]
            cx = self.allLines[index][20:28] 
            cy = self.allLines[index][28:36]
            cz = self.allLines[index][36:44]
            
            x = float(cx)*10
            y = float(cy)*10
            z = float(cz)*10
            atom._coords = [[x, y, z]]
             
            atom._charges = []
            atom.segID =  mol.name   
            atom.normalname = name
            atom.number = int(self.allLines[index][15:20])
            atom.elementType = name[0]
            mol.atmNum[atom.number] = atom
            atom.altname = None    
            atom.hetatm = 0
        mol.name = os.path.split(os.path.splitext(self.filename)[0])[-1]
        mol.allAtoms = mol.chains.residues.atoms
        mol.parser = self
        mol.levels = [Protein, Chain, Residue, Atom]
        name = ''
        for n in molList.name:
            name = n + ','
        name = name[:-1]
        molList.setStringRepr( name )
        strRpr = name + ':::'
        molList.allAtoms.setStringRepr( strRpr )
        for m in molList:
            mname = m.name
            strRpr = mname + ':::'
            m.allAtoms.setStringRepr( strRpr )
            strRpr = mname + ':'
            m.chains.setStringRepr( strRpr )
            for c in m.chains:
                cname = c.id
                strRpr = mname + ':' + cname + ':'
                c.residues.setStringRepr( strRpr )
                for r in c.residues:
                    rname = r.name
                    strRpr = mname + ':' + cname + ':' + rname + ':'
                    r.atoms.setStringRepr( strRpr )        
        return molList
        

    def getMoleculeInformation(self):
        """ Function to retrieve the general informations on the molecule.
        This information is used by the molecule chooser to provide
        informations on the molecule selected.
        """
        molStr = ''
        return molStr
    
    def configureProgressBar( self, **kw ):
        # this method is to be implemented by the user from outside
        pass
        
if __name__ == '__main__':
    parser = groParser( filename='/usr/share/gromacs/tutor/methanol/conf.gro' )
    print "Reading molecule"
    mol = parser.parse()
