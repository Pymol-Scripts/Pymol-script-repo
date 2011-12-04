########################################################################
#
# Date: April 2006 Authors: Guillaume Vareille, Michel Sanner
#
#    sanner@scripps.edu
#    vareille@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
# revision:
#
#########################################################################
#
# $Header$
#
# $Id$
#

from NetworkEditor.datatypes import AnyArrayType

class MSMSobjectType(AnyArrayType):

    from mslib import MSMS
    def __init__(self, name='MSMSobject', datashape=None, color='#6060c3',
                 shape='pentagon', width=None, height=None, 
                 dataDescr='MSMS Instance',
                 klass=MSMS):
        
        AnyArrayType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)



class CrystToReal(AnyArrayType):

    from mglutil.math.crystal import Crystal
    def __init__(self, name='CrystToReal', datashape=None, color='#fe9200',
                 shape='oval', width=None, height=None, 
                 dataDescr='Crystal to Real conversion object',
                 klass=Crystal):
        
        AnyArrayType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)



class TreeNodeType(AnyArrayType):

    from MolKit.tree import TreeNode
    def __init__(self, name='TreeNode', datashape=None, color='#c64e70',
                 shape='oval', width=8, height=12, dataDescr='TreeNode Instance',
                 klass=TreeNode):
        
        AnyArrayType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)



class ChainType(TreeNodeType):

    def __init__(self):
        TreeNodeType.__init__(self)
        self.data['name'] = 'Chain'
        self.data['dataDescr'] = 'Chain Instance'
        self.data['color'] = '#e64e70'
        self.data['shape'] = 'oval'

    from MolKit.protein import Chain
    def __init__(self, name='Chain', datashape=None, color='#e64e70',
                 shape='oval', width=8, height=12, dataDescr='Chain Instance',
                 klass=Chain):
        
        TreeNodeType.__init__(self, name=name, datashape=datashape, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, dataDescr=dataDescr)



class AtomType(TreeNodeType):

    from MolKit.molecule import Atom
    def __init__(self, name='Atom', datashape=None, color='#fe92a0',
                 shape='oval', width=8, height=12, dataDescr='Atom Instance',
                 klass=Atom):
        
        TreeNodeType.__init__(self, name=name, datashape=datashape, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, dataDescr=dataDescr)



class MoleculeType(TreeNodeType):

    from MolKit.molecule import Molecule
    def __init__(self, name='Molecule', datashape=None, color='#c64e70',
                 shape='oval', width=8, height=12, dataDescr='Molecule Instance',
                 klass=Molecule):
        
        TreeNodeType.__init__(self, name=name, datashape=datashape, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, dataDescr=dataDescr)



class ResidueType(TreeNodeType):

    from MolKit.protein import Residue
    def __init__(self, name='Residue', datashape=None, color='#f64e70',
                 shape='oval', width=8, height=12, dataDescr='Residue Instance',
                 klass=Residue):
        
        TreeNodeType.__init__(self, name=name, datashape=datashape, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, dataDescr=dataDescr)



class TreeNodeSetType(AnyArrayType):

    from MolKit.tree import TreeNodeSet
    def __init__(self, name='TreeNodeSet', datashape=None, color='#c64e70',
                 shape='oval', width=None, height=None, dataDescr='TreeNodeSet Instance',
                 klass=TreeNodeSet):
        
        AnyArrayType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)


    def cast(self, data):
        """returns a success status (true, false) and the coerced data
"""
        from MolKit.tree import TreeNode
        if self['datashape'] is None: 
            if isinstance(data, TreeNode):
                return True, data.setClass([data])
        else:
            lArray = Numeric.array(data)
            lArray0 = lArray[0]
            while hasattr(lArray0,'shape'):
                lArray0 = lArray0[0]
            if isinstance(lArray0, TreeNode):
                # this is only casting the first elt of the array (its better then nothing)
                return True, lArray0.setClass([lArray0])
        return False, None





class MoleculeSetType(TreeNodeSetType):

    from MolKit.molecule import MoleculeSet
    def __init__(self, name='MoleculeSet', datashape=None, color='#c64e70',
                 shape='oval', width=None, height=None, dataDescr='MoleculeSet Instance',
                 klass=MoleculeSet):
        
        TreeNodeSetType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)




class ChainSetType(TreeNodeSetType):

    from MolKit.protein import ChainSet
    def __init__(self, name='ChainSet', datashape=None, color='#e64e70',
                 shape='oval', width=None, height=None, dataDescr='ChainSet Instance',
                 klass=ChainSet):
        
        TreeNodeSetType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)




class ResidueSetType(TreeNodeSetType):

    from MolKit.protein import ResidueSet
    def __init__(self, name='ResidueSet', datashape=None, color='#f64e70',
                 shape='oval', width=None, height=None, dataDescr='ResidueSet Instance',
                 klass=ResidueSet):
        
        TreeNodeSetType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)




class AtomSetType(TreeNodeSetType):

    from MolKit.molecule import AtomSet
    def __init__(self, name='AtomSet', datashape=None, color='#fe92a0',
                 shape='oval', width=None, height=None, dataDescr='AtomSet Instance',
                 klass=AtomSet):
        
        TreeNodeSetType.__init__(self, name=name, color=color, 
                              shape=shape, width=width, height=height, 
                              klass=klass, datashape=datashape, 
                              dataDescr=dataDescr)

