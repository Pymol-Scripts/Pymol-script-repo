#########################################################################
#
# Date: Apr 2006 Authors: Ruth Huey, Garrett Morris
#
#    rhuey@scripps.edu
#    garrett@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Ruth Huey, Garrett Morris and TSRI
#
#########################################################################

import warnings

from Vision import UserLibBuild
from NetworkEditor.items import NetworkNode
from MolKit.molecule import Atom, AtomSet, Molecule, MoleculeSet
from MolKit.protein import Residue, ResidueSet, Chain, ChainSet
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from AutoDockTools.atomTypeTools import NonpolarHydrogenMerger, LonepairMerger
from AutoDockTools.atomTypeTools import AromaticCarbonManager, SolvationParameterizer
from AutoDockTools.MoleculePreparation import RotatableBondManager
from AutoDockTools.MoleculePreparation import LigandWriter, AD4LigandWriter



def importAdtLib(net):
    try:
        from AutoDockTools.VisionInterface.AdtNodes import adtlib
        net.editor.addLibraryInstance(
            adtlib, 'AutoDockTools.VisionInterface.AdtNodes', 'adtlib')
    except:
        warnings.warn(
            'Warning! Could not import adtlib from AutoDockTools.VisionInterface')
            

class GridParameterFileBrowserNE(NetworkNode):
    """A specialized Tkinter Filebrowser. Double-clicking into the entry opens the
filebrowser."""
    
    def __init__(self, name='Grid Parameter File Browser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        #self.readOnly = 1
        code = """def doit(self, filename):
    if filename:
        self.outputData(filename=filename)
"""

        self.setFunction(code)

        # show the entry widget by default
        self.inNodeWidgetVisibleByDefault = True

        fileTypes=[('gpf', '*')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'read file', 'width':16,
            'labelCfg':{'text':'gpf file:'},
            }


        #self.widgetDescr['filename'] = {
        #    'class':'NEEntryWithFileBrowser', 'master':'node', 'width':16,
        #    'initialValue':'', 'lockedOnPort':True, 
        #    'labelCfg':{'text':'Filename: '}
        #    }

        self.inputPortsDescr.append(datatype='string', name='filename')

        self.outputPortsDescr.append(datatype='string', name='filename')


class DockingParameterFileBrowserNE(NetworkNode):
    """A specialized Tkinter Filebrowser. Double-clicking into the entry opens the
filebrowser."""
    
    def __init__(self, name='Docking Parameter File Browser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        #self.readOnly = 1
        code = """def doit(self, filename):
    if filename:
        self.outputData(filename=filename)
"""

        self.setFunction(code)

        # show the entry widget by default
        self.inNodeWidgetVisibleByDefault = True

        fileTypes=[('dpf', '*')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'read file', 'width':16,
            'labelCfg':{'text':'dpf file:'},
            }

        self.inputPortsDescr.append(datatype='string', name='filename')

        self.outputPortsDescr.append(datatype='string', name='filename')



class ReadGridParameterFile(NetworkNode):
    #mRequiredTypes = {}
    #mRequiredSynonyms = [
    #]
    def __init__(self, constrkw = {},  name='ReadGridParameterFile', **kw):
        kw['constrkw'] = constrkw
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)
        fileTypes=[('gpf', '*')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'read file', 'width':16,
            'labelCfg':{'text':'file:'},
            }

        code = """def doit(self, template_gpf_filename):
    if template_gpf_filename:
            from AutoDockTools.GridParameters import GridParameters
            gpo = GridParameters()
            gpo.read(template_gpf_filename)
            self.outputData(gpo=gpo)
"""
        self.configure(function=code)

        self.inputPortsDescr.append(
            {'name': 'template_gpf_filename', 'cast': True, 'datatype': 'string', 'balloon': 'template grid parameter filename', 'required': False, 'height': 8, 'width': 12, 'shape': 'oval', 'color': 'white'})
        self.outputPortsDescr.append(
            {'name': 'gpo', 'datatype': 'None', 'balloon': 'gpo,  grid parameter object,  instance of AutoDockTools.GridParameters', 'height': 8, 'width': 12, 'shape': 'diamond', 'color': 'white'})
        self.widgetDescr['template_gpf_filename'] = {
            'initialValue': '', 'labelGridCfg': {'column': 0, 'row': 0}, 'master': 'node', 'widgetGridCfg': {'labelSide': 'left', 'column': 1, 'row': 0}, 'labelCfg': {'text': ''}, 'class': 'NEEntryWithFileBrowser'}

    def beforeAddingToNetwork(self, net):
        try:
            ed = net.getEditor()
        except:
            import traceback; traceback.print_exc()
            print 'Warning! Could not import widgets'



class ReadDockingParameterFile(NetworkNode):
    #mRequiredTypes = {}
    #mRequiredSynonyms = [
    #]
    def __init__(self, constrkw = {},  name='ReadDockingParameterFile', **kw):
        kw['constrkw'] = constrkw
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)
        fileTypes=[('dpf', '*')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'read file', 'width':16,
            'labelCfg':{'text':'file:'},
            }

        code = """def doit(self, template_dpf_filename):
    if template_dpf_filename:
            from AutoDockTools.DockingParameters import DockingParameters
            gpo = DockingParameters()
            gpo.read(template_dpf_filename)
            self.outputData(gpo=gpo)
"""
        self.configure(function=code)

        self.inputPortsDescr.append(
            {'name': 'template_dpf_filename', 'cast': True, 'datatype': 'string', 'balloon': 'template grid parameter filename', 'required': False, 'height': 8, 'width': 12, 'shape': 'oval', 'color': 'white'})
        self.outputPortsDescr.append(
            {'name': 'gpo', 'datatype': 'None', 'balloon': 'gpo,  grid parameter object,  instance of AutoDockTools.DockingParameters', 'height': 8, 'width': 12, 'shape': 'diamond', 'color': 'white'})
        self.widgetDescr['template_dpf_filename'] = {
            'initialValue': '', 'labelGridCfg': {'column': 0, 'row': 0}, 'master': 'node', 'widgetGridCfg': {'labelSide': 'left', 'column': 1, 'row': 0}, 'labelCfg': {'text': ''}, 'class': 'NEEntryWithFileBrowser'}

    def beforeAddingToNetwork(self, net):
        try:
            ed = net.getEditor()
        except:
            import traceback; traceback.print_exc()
            print 'Warning! Could not import widgets'


###class ReadGridParameterFile(NetworkNode):
###    """Read a Grid Parameter file [using Python's readlines() command.]
###Double-clicking on the node opens a text entry widget to type the file name.
###In addition, double-clicking in the text entry opens a file browser window.

###Input Ports
###    filename: name of the file to be opened

###Output Ports
###    data: a list of strings
###"""

###    def __init__(self, name='Read Parameter File', **kw):
###        kw['name'] = name
###        apply( NetworkNode.__init__, (self,), kw )

###        #self.readOnly = 1
###        code = """def doit(self, filename):
###    if filename and len(filename):
###        gpo = GridParameters()
###        gpo.read(filename)
###        #f = open(filename)
###        #datastream = f.readlines()
###        #f.close()
###        #if datastream:
###        self.outputData(data=gpo)
###"""

###        self.setFunction(code)

###        fileTypes=[('gpf', 'gpf')]

###        self.widgetDescr['filename'] = {
###            'class':'NEEntryWithFileBrowser', 'master':'node',
###            'filetypes':fileTypes, 'title':'read file', 'width':16,
###            'labelCfg':{'text':'file:'},
###            }


###        self.inputPortsDescr.append(datatype='string', name='filename')
###        self.outputPortsDescr.append(datatype='instance', name='gpo')


class RemoveWaters(NetworkNode):
    """ removes water residues
Process entails:
#    1.  looping over each molecule in input
#    2.  selecting all atoms in water residues
#    3.  removing all bonds from each atom
#    4.  removing each atom from its parent residue
#    5.  removing parent residue if it has no remaining atoms
#    6.  removing parent chain if it has no remaining residues
#    7.  resetting allAtoms attribute of molecule
#    8.  returning number of water residues removed
Input:  molecules (MoleculeSet)
Output: molecules where all atoms in HOH residues have been removed(MoleculeSet)"""
    
    def __init__(self, name='RemoveWaters', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')
        ip.append(datatype='str', required=False, name='residue_type_str', defaultValue='HOH')
        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='molecules_with_no_water_residues')
        op.append(datatype='int', name='num_water_res')

        code = """def doit(self, molecules, residue_type_str):
    if molecules:
        from MolKit.molecule import BondSet
        lenHOHs = 0
        for mol in molecules:
            hohs = mol.allAtoms.get(lambda x: x.parent.type==residue_type_str)
            if hohs:
                #remove(hohs)
                lenHOHs = len(hohs)
                for h in hohs:
                    for b in h.bonds:
                        c = b.atom1
                        if c==h:
                            c = b.atom2
                        c.bonds.remove(b)
                    h.bonds = BondSet()
                    res = h.parent
                    h.parent.remove(h)
                    if len(h.parent.children)==0:
                        res = h.parent
                        chain = res.parent
                        #print 'removing residue: ', res.name
                        chain.remove(res)
                        if len(chain.children)==0:
                            mol = chain.parent
                            print 'removing chain', chain.id
                            mol.remove(chain)
                            del chain
                        del res
                    del h 
                #fix allAtoms short cut
                mol.allAtoms = mol.chains.residues.atoms    
    self.outputData(molecules_with_no_water_residues=molecules, num_water_res=lenHOHs)"""
        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass




class MergeNonPolarHydrogens(NetworkNode):
    """ merges nonpolar hydrogens 
Process entails:
#    1.  adding charge on nonpolar hydrogen to charge of carbon atom to which it is bonded
#    2.  removing the nonpolar hydrogen from the molecule
#    3.  resetting allAtoms attribute of molecule
#    4.  returning number of nonpolar hydrogens removed
Input:  mols (MoleculeSet)
Output: mols where each non-polar hydrogen atom has been removed(MoleculeSet)"""
    
    def __init__(self, name='Nonpolar Hydrogen Merger', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')
        ip.append(datatype='int', required=False, name='renumber', defaultValue=1)

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='mols')
        op.append(datatype='int', name='num_nphs')

        code = """def doit(self,mols, renumber):
    if mols:
        nph_merger = NonpolarHydrogenMerger()
        num_nphs = 0
        for mol in mols:
            if not len(mol.allAtoms.bonds[0]):
                mol.buildBondsByDistance()
            num_nphs = nph_merger.mergeNPHS(mol.allAtoms, renumber=renumber)
        self.outputData(mols=mols, num_nphs=num_nphs)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class  MergeLonePairs(NetworkNode):
    """ merges lone pairs
Process entails:
#    1.  adding charge on lone pair 'atom' to charge of carbon atom to which it is 'bonded'
#    2.  removing the lone pair 'atom' from the molecule
#    3.  resetting allAtoms attribute of molecule
#    4.  returning number of lone pair 'atoms' removed
Input:  mols (MoleculeSet)
Output: mols where each lone pair 'atom' has been removed(MoleculeSet)"""
    
    def __init__(self, name='Lone Pair Merger', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')
        ip.append(datatype='int', required=False, name='renumber', defaultValue=1)

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='mols')
        op.append(datatype='int', name='num_lps')

        code = """def doit(self,mols, renumber):
    if mols:
        lps_merger = LonepairMerger()
        num_lps = 0
        for mol in mols:
            if not len(mol.allAtoms.bonds[0]):
                mol.buildBondsByDistance()
            lps = lps_merger.mergeLPS(mol.allAtoms, renumber=renumber)
            num_lps += len(lps)
        self.outputData(mols=mols, num_lps=num_lps)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass


class  ManageAromaticCarbons(NetworkNode):
    """ manages assigning autodock carbon atom types: aliphatic and aromatic
Process entails:
#    1.  'setAromaticCarbons' method detects cyclic aromatic carbons: renaming 
#    them 'A' and setting autodock_element to 'A'. Cyclic carbons are 'aromatic' 
#    if the angle between adjacent normals  is  less than a specified 'cutoff' 
#    angle which is 7.5 degrees by default. Returns atoms which are changed 
#    2.  provides widget for changing the 'cutoff'
#    3.  NOT_IMPLEMENTED:provides method 'set_carbon_names' for forcing any 'C' to 'A' and
#    the opposite. Returns atoms which are changed
Input:  mols (MoleculeSet)
Output: mols where aromatic carbons have been detected (MoleculeSet)"""
    
    def __init__(self, name='Manage Aromatic Carbons', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['cutoff'] = {
                'class':'NEDial', 'size':50,
                'oneTurn':5.0, 'min':1.0, 'lockMin':1, 'type':'float',
                'initialValue':7.5,
                'labelGridCfg':{'sticky':'w'},
                'labelCfg':{'text':'cutoff'},
                }

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')
        ip.append(datatype='float', required=False, name='cutoff', defaultValue=7.5)

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='mols')
        op.append(datatype='int', name='num_aromaticCs')

        code = """def doit(self, mols, cutoff):
    if mols:
        aromC_manager = AromaticCarbonManager(cutoff=cutoff)
        num_aromCs = 0
        for mol in mols:
            if not len(mol.allAtoms.bonds[0]):
                mol.buildBondsByDistance()
            aromCs = aromC_manager.setAromaticCarbons(mol, cutoff=cutoff)
            num_aromCs+=len(aromCs)
        self.outputData(mols=mols, num_aromaticCs=num_aromCs)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass


class  Assign_AD4Types(NetworkNode):
    """ assigns autodock4-style 'autodock_element' to atoms
Process entails:
#    1.  distinguishing between nitrogens which do accept hydrogen-bonds and
#    those which do not (because they already have bonds to hydrogens)
#    2.  distinguishing between oxygens which do accept hydrogen-bonds and
#    those which do not (because they already have bonds to hydrogens)
#    3.  setting autodock_element to 'A' for carbon atoms in cycles in standard amino acids
#    4.  setting autodock_element to 'HD' for all hydrogen atoms   
#    5.  setting autodock_element for all other atoms to the atom's element
#    NOTE: more types can be added if more distinctions are supported by
#    autodock
Input:  mols (MoleculeSet)
Output: typed_mols where each atom has autodock_element field(MoleculeSet)"""
    
    def __init__(self, name='AD4_typer', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')
        ip.append(datatype='int', required=False, name='set_aromatic_carbons', defaultValue=1)
        ip.append(datatype='int', required=False, name='reassign', defaultValue=1)
        ip.append(datatype='int', required=False, name='renameAtoms', defaultValue=0)

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='typed_mols')

        code = """def doit(self,mols, set_aromatic_carbons=1, reassign=1, renameAtoms=0):
    if mols:
        at_typer = AutoDock4_AtomTyper(set_aromatic_carbons=1, renameAtoms=renameAtoms)
        for mol in mols:
            if not len(mol.allAtoms.bonds[0]):
                mol.buildBondsByDistance()
            at_typer.setAutoDockElements(mol, reassign=reassign)
        self.outputData(typed_mols=mols)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class  Add_SolvationParameters(NetworkNode):
    """ assigns autodock3-style 'solvation parameters' to atoms
Process entails:
#    1.  distinguishing between aromatic and aliphatic carbons
#    2.  look-up table of solvation volumes (AtSolVol)
#    3.  setting AtSolPar and AtVol
Input:  mols (MoleculeSet)
Output: typed_mols where each atom has SolVol and AtSolPar fields(MoleculeSet)"""
    
    def __init__(self, name='AD4_SolvationParameterizer', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='typed_mols')
        op.append(datatype='list', name='AtSolPar')

        code = """def doit(self,mols):
    if mols:
        SP = SolvationParameterizer()
        for mol in mols:
            unknown_atoms = SP.addParameters(mol.chains.residues.atoms)
            #?keep information about unknown_atoms?
            if unknown_atoms is not None: mol.unknown_atoms = len(unknown_atoms)
        self.outputData(typed_mols=mols, AtSolPar=mols.allAtoms.AtSolPar)\n"""
        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class  ManageRotatableBonds(NetworkNode):
    """ manages setting flexibility pattern in molecule
Process entails:
#    1.  distinguishing between possibly-rotatable and non-rotatable bonds
#    2.  turning on/off classes of possibly-rotatable bonds such as 
#           amide, guanidinium, peptide-backbone
#    3.  optionally turning on/off specific bonds between named atoms
#    4.  optionally limiting the total number of rotatable bonds 
#           toggling either: 
#               those which move the most atoms
#                   or 
#               those which move the fewest atoms
Input:  mol (Molecule)
Output: mol with marked bonds """
    
    def __init__(self, name='Manage Rotatable Bonds', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols', defaultValue='auto')
        ip.append(datatype='string', required=False, name='root')
        ip.append(datatype='string', required=False, name='allowed_bonds', defaultValue='backbone')
        ip.append(datatype='int', required=False, name='check_for_fragments', defaultValue=0)
        ip.append(datatype='string', required=False, name='bonds_to_inactivate', defaultValue='')
        ip.append(datatype='string', required=False, name='limit_torsions', defaultValue='')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='mols')

        code = """def doit(self, mols, root, allowed_bonds, check_for_fragments,
bonds_to_inactivate, limit_torsions):
    if mols:
        #mol = mols[0]
        for mol in mols:
            if not len(mol.allAtoms.bonds[0]):
                mol.buildBondsByDistance()
            print "root=", root
            mol.RBM = RotatableBondManager(mol, allowed_bonds, root,
                            check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate)
            if limit_torsions:
                mol.RBM.limit_torsions(limit_torsions)
        self.outputData(mols=mols)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class  Ligand_Writer(NetworkNode):
    """ writes autodock3 ligand pdbq file
Process entails:
#    1.  writing REMARK records about torsions
#    2.  writing ROOT/ENDROOT records about rigid portion of ligand
#    3.  writing BRANCH/ENDBRANCH records about movable sections of ligand
#    4.  writing TORSDOF record showing torsional degrees of freedom
Input:  mol (Molecule), output_filename (string)
Output: mol """
    
    def __init__(self, name='Ligand_Writer', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        fileTypes=[('pdbq', '*.pdbq'), ('all', '*')]

        self.widgetDescr['output_filename'] = {
            'class':'NEEntryWithFileSaver', 'master':'node',
            'filetypes':fileTypes, 'title':'save AD3 Ligand', 'width':10,
            'labelCfg':{'text':'file:'},
            }
        
        ip = self.inputPortsDescr
        ip.append(datatype='Molecule', name='mol')
        ip.append(datatype='string',  name='output_filename')

        op = self.outputPortsDescr
        op.append(datatype='Molecule', name='mol')

        code = """def doit(self, mol, output_filename):
    if mol:
        #mol = mols[0]
        #check for bonds with 'possibleTors/activeTOrs keys'
        #check for root
        #check for TORSDOF
        #check for charges
        writer = LigandWriter()
        writer.write(mol, output_filename)
        self.outputData(mol=mol)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class  Ligand_Writer_AD4(NetworkNode):
    """ writes autodock4 ligand pdbqt file
Process entails:
#    1.  writing REMARK records about torsions
#    2.  writing ROOT/ENDROOT records about rigid portion of ligand
#    3.  writing BRANCH/ENDBRANCH records about movable sections of ligand
#    4.  writing TORSDOF record showing torsional degrees of freedom
Input:  mol (Molecule), output_filename (string)
Output: mol, output_filename """
    
    def __init__(self, name='Ligand_Writer_AD4', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        fileTypes=[('pdbqt', '*.pdbqt'), ('all', '*')]

        self.widgetDescr['output_filename'] = {
            'class':'NEEntryWithFileSaver', 'master':'node',
            'filetypes':fileTypes, 'title':'save AD4 Ligand', 'width':10,
            'labelCfg':{'text':'file:'},
            }
        
        ip = self.inputPortsDescr
        ip.append(datatype='Molecule', name='mol')
        ip.append(datatype='string', name='output_filename')

        op = self.outputPortsDescr
        op.append(datatype='Molecule', name='mol')
        op.append(datatype='string', name='output_filename')

        code = """def doit(self, mol, output_filename):
    if mol:
        writer = AD4LigandWriter()
        writer.write(mol, output_filename)
        self.outputData(mol=mol, output_filename=output_filename)\n"""
        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



####class AdtPrepareLigand(NetworkNode):
####    """ formats ligand for autodock3 using AutoDockTools.MoleculePreparation.LigandPreparation.
####Process entails:
####    1. possible clean_up:
####            removing lone pairs
####            merging non_polar hydrogens
####            adding bonds to atoms with no bonds
####    2. making atoms in ligand conform to autodock3 atom types:
####            distinction between carbons and cyclic-aromatic carbons, 
####            no non-polar hydrogens
####    3. adding partial charges (gasteiger by default)
####    4. establishing 'torTree' for flexibility pattern by setting 'root' and 'rotatable' bonds
####    5. writing 'pdbq' file
####Input:  mols (MoleculeSet)
####Output: AD3ligand (Molecule)"""


###class Adt4PrepareLigand(NetworkNode):
###    """ formats ligand for autodock4 using AutoDockTools.MoleculePreparation.AD4LigandPreparation.
###Process entails:
###    1. possible clean_up:
###            removing lone pairs
###            merging non_polar hydrogens
###            adding bonds to atoms with no bonds
###    2. making atoms in ligand conform to autodock3 atom types:
###            distinction between carbons and cyclic-aromatic carbons, 
###            no non-polar hydrogens
###    3. adding partial charges (gasteiger by default)
###    4. establishing 'torTree' for flexibility pattern by setting 'root' and 'rotatable' bonds
###    5. writing 'pdbqt' file
###Input:  mols (MoleculeSet)
###Output: AD4ligand (Molecule)"""


###class AdtPrepareReceptor(NetworkNode):
###    """ formats receptor for autodock3 using AutoDockTools.MoleculePreparation.ReceptorPreparation.
###Process entails:
###    1. possible clean_up:
###            removing lone pairs
###            merging non_polar hydrogens
###            adding bonds to atoms with no bonds
###    2. making atoms in receptor conform to autodock3 atom types:
###            distinction between carbons and cyclic-aromatic carbons, 
###            polar hydrogens but no non-polar hydrogens
###    3. adding partial charges (Kollman by default)
###    4. writing 'pdbqs' file
###Input:  mols (MoleculeSet)
###Output: AD3receptor (Molecule)"""

###class Adt4PrepareReceptor(NetworkNode):
###    """ formats receptor for autodock4 using AutoDockTools.MoleculePreparation.AD4ReceptorPreparation.
###Process entails:
###    1. possible clean_up:
###            removing lone pairs
###            merging non_polar hydrogens
###            adding bonds to atoms with no bonds
###    2. making atoms in receptor conform to autodock4 atom types:
###            distinction between carbons and cyclic-aromatic carbons, 
###            distinction between hydrogen-bonding and non-hydrogen-bonding nitrogens,
###            no non-polar hydrogens
###    3. adding partial charges (gasteiger by default)
###    4. writing 'pdbqt' file
###Input:  mols (MoleculeSet)
###Output: AD4receptor (Molecule)"""
    

class AdtPrepareGpf3(NetworkNode):
    """ writes parameter file for autogrid3
Process entails:
    1.  setting ligand
    2.  setting receptor
    3.  setting specified values of various parameters
    4.  writing gpf file for autogrid3
Input:  ligand_filename, receptor_filename, optional parameter dictionary
Output: gpf ('string')"""
    
    def __init__(self, name='Prepare Autogrid3 Gpf', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='ligand_filename')
        ip.append(datatype='string', name='receptor_filename')
        ip.append(datatype='string', required=False, name='gpf_filename', defaultValue='')
        ip.append(datatype='dict', required=False, name='parameters', defaultValue={})
        ip.append(datatype='string', required=False, name='outputfilename', defaultValue='')

        op = self.outputPortsDescr
        op.append(datatype='string', name='ag3_parameter_file')

        code = """def doit(self, ligand_filename, receptor_filename, gpf_filename, parameters, outputfilename):
    if ligand_filename and receptor_filename:
        from AutoDockTools.GridParameters import GridParameterFileMaker
        gpfm = GridParameterFileMaker()
        gpfm.set_ligand(ligand_filename)
        gpfm.set_receptor(receptor_filename)
        if gpf_filename:
            gpfm.read_reference(gpf_filename)
        if len(parameters):
            gpfm.set_grid_parameters(parameters)
        if not outputfilename:
            outputfilename = gpfm.ligand.name+'_'+gpfm.receptor_stem + ".gpf"
        gpfm.write_gpf(outputfilename)
        self.outputData(ag3_parameter_file=outputfilename)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class AdtPrepareGpf4(NetworkNode):
    """ writes parameter file for autogrid4
Process entails:
    1.  setting ligand
    2.  setting receptor
    3.  setting specified values of various parameters
    4.  writing gpf file for autogrid4
Input:  ligand_filename, receptor_filename, optional parameter dictionary
Output: gpf ('string')"""
    
    def __init__(self, name='Prepare Autogrid4 Gpf', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='ligand_filename')
        ip.append(datatype='string', name='receptor_filename')
        ip.append(datatype='string', required=False, name='gpf_filename', defaultValue='')
        ip.append(datatype='dict', required=False, name='parameters', defaultValue={})
        ip.append(datatype='string', required=False, name='outputfilename', defaultValue='')

        op = self.outputPortsDescr
        op.append(datatype='string', name='ag4_parameter_file')

        code = """def doit(self, ligand_filename, receptor_filename, gpf_filename, parameters, outputfilename):
    if ligand_filename and receptor_filename:
        from AutoDockTools.GridParameters import GridParameter4FileMaker
        gpfm = GridParameter4FileMaker()
        gpfm.set_ligand(ligand_filename)
        gpfm.set_receptor(receptor_filename)
        if gpf_filename:
            gpfm.read_reference(gpf_filename)
        if len(parameters):
            gpfm.set_grid_parameters(parameters)
        if not outputfilename:
            outputfilename = gpfm.ligand.name+'_'+gpfm.receptor_stem + ".gpf"
        gpfm.write_gpf(outputfilename)
        self.outputData(ag4_parameter_file=outputfilename)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class AdtPrepareDpf3(NetworkNode):
    """ writes parameter file for autodock3
Process entails:
    1.  setting ligand
    2.  setting receptor
    3.  setting specified values of various parameters
    4.  writing dpf file for autogrid3
Input:  ligand_filename, receptor_filename, optional parameter dictionary
Output: dpf ('string')"""
    
    def __init__(self, name='Prepare Autodock3 Dpf', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='ligand_filename')
        ip.append(datatype='string', name='receptor_filename')
        ip.append(datatype='string', required=False, name='dpf_filename', defaultValue='')
        ip.append(datatype='dict', required=False, name='parameters', defaultValue={})
        ip.append(datatype='string', required=False, name='outputfilename', defaultValue='')

        op = self.outputPortsDescr
        op.append(datatype='string', name='ad3_parameter_file')

        code = """def doit(self, ligand_filename, receptor_filename, dpf_filename, parameters, outputfilename):
    if ligand_filename and receptor_filename:
        from AutoDockTools.DockingParameters import DockingParameterFileMaker
        dpfm = DockingParameterFileMaker()
        dpfm.set_ligand(ligand_filename)
        dpfm.set_receptor(receptor_filename)
        if dpf_filename:
            dpfm.read_reference(dpf_filename)
        if len(parameters):
            dpfm.set_docking_parameters(parameters)
        if not outputfilename:
            outputfilename = dpfm.ligand.name+'_'+dpfm.receptor_stem + ".dpf"
        dpfm.write_dpf(outputfilename)
        self.outputData(ad3_parameter_file=outputfilename)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass



class AdtPrepareDpf4(NetworkNode):
    """ writes parameter file for autodock4
Process entails:
    1.  setting ligand
    2.  setting receptor
    3.  setting specified values of various parameters
    4.  writing dpf file for autodock4
Input:  ligand_filename, receptor_filename, optional parameter dictionary
Output: dpf ('string')"""
    
    def __init__(self, name='Prepare Autodock4 Dpf', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='ligand_filename')
        ip.append(datatype='string', name='receptor_filename')
        ip.append(datatype='string', required=False, name='dpf_filename', defaultValue='')
        ip.append(datatype='dict', required=False, name='parameters', defaultValue={})
        ip.append(datatype='string', required=False, name='outputfilename', defaultValue='')

        op = self.outputPortsDescr
        op.append(datatype='string', name='ad4_parameter_file')

        code = """def doit(self, ligand_filename, receptor_filename, dpf_filename, parameters, outputfilename):
    if ligand_filename and receptor_filename:
        from AutoDockTools.DockingParameters import DockingParameter4FileMaker
        dpfm = DockingParameter4FileMaker()
        dpfm.set_ligand4(ligand_filename)
        dpfm.set_receptor4(receptor_filename)
        if dpf_filename:
            dpfm.read_reference4(dpf_filename)
        if len(parameters):
            dpfm.set_docking_parameters(parameters)
        if not outputfilename:
            outputfilename = dpfm.ligand.name+'_'+dpfm.receptor_stem + ".dpf"
        dpfm.write_dpf(outputfilename)
        self.outputData(ad4_parameter_file=outputfilename)\n"""

        self.setFunction(code)


    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass


#class AdtPrepareFlexDocking4
#class AdtPdbqToPdbqt
#class AdtPdbqsToPdbqt
#class AdtGpf3ToGpf4
#class AdtDpf3ToDpf4
#class AdtRunAutogrid3
#class AdtRunAutogrid4
#class AdtRunAutodock3
#class AdtRunAutodock4
#class AdtSummarizeResults
#class AdtSummarizeResults4
#class AdtSummarizeXMLResults4
#class AdtPixelMapResults
#objects:
#MolecularPreparation classes
#   AutoDockBondClassifier
#   ReceptorWriter<==MolKitNodes/WriteMolecule
#   AD4ReceptorWriter<==MolKitNodes/WriteMolecule
#AD4FlexibleDockingPreparation
#prepare_ligand_dict
#ReceptorWriter
#AD4ReceptorWriter
#LigandPreparation
#AD4LigandPreparation
#AutoDock3_AtomTyper
#DockingParameters
#DockingParameterFileMaker
#GridParameters
#GridParameterFileMaker
#??Docking??


from Vision.VPE import NodeLibrary
adtlib = NodeLibrary('adtlib', '#4444FF')

#macros from other files in this directory
#from AutoDockTools.VisionInterface.PrepareAD4Molecule import PrepareAD4Molecule
#adtlib.addNode(PrepareAD4Molecule, 'Prepare AD4Molecule', 'Macro')
#from AutoDockTools.VisionInterface.AD4Ligands import AD4Ligands
#adtlib.addNode(AD4Ligands, 'Prepare AD4 Ligands', 'Macro')
#from AutoDockTools.VisionInterface.Prepare_AD4_Ligands import Prepare_AD4_Ligands
#adtlib.addNode(Prepare_AD4_Ligands, 'Prepare AD4 Ligands', 'Macro')
###from AutoDockTools.VisionInterface.PrepareAD3Ligand import PrepareAD3Ligand
###adtlib.addNode(PrepareAD3Ligand, 'Prepare AD3 Ligand', 'Macro')
###from AutoDockTools.VisionInterface.PrepareAD3Receptor import PrepareAD3Receptor
###adtlib.addNode(PrepareAD3Receptor, 'Prepare AD3Receptor', 'Macro')
###from AutoDockTools.VisionInterface.AD3Gpf import AD3Gpf
###adtlib.addNode(AD3Gpf, 'Prepare AD3 Gpf ', 'Macro')
###from AutoDockTools.VisionInterface.AD3Dpf import AD3Dpf
###adtlib.addNode(AD3Dpf, 'Prepare AD3 Dpf ', 'Macro')
#from AutoDockTools.VisionInterface.GPF4 import GPF4
#adtlib.addNode(GPF4, 'Prepare AD4 Gpf ', 'Macro')
from AutoDockTools.VisionInterface.Docking import Docking
adtlib.addNode(Docking, 'Docking', 'Macro')
from AutoDockTools.VisionInterface.recluster import recluster
adtlib.addNode(recluster, 'recluster...', 'Macro')

adtlib.addNode(GridParameterFileBrowserNE, 'Grid Parameter File Browser', 'Input')
adtlib.addNode(DockingParameterFileBrowserNE, 'Docking Parameter File Browser', 'Input')

#adtlib.addNode(ReadGridParameterFile, 'Read Grid Parameter File', 'Input')
#adtlib.addNode(ReadDockingParameterFile, 'Read Docking Parameter File', 'Input')


#adtlib.addNode(AdtPrepareGpf3, 'Prepare AD3Gpf', 'Macro')
#adtlib.addNode(AdtPrepareGpf4, 'Prepare AD4Gpf', 'Macro')
#adtlib.addNode(AdtPrepareDpf3, 'Prepare AD3Dpf', 'Macro')
#adtlib.addNode(AdtPrepareDpf4, 'Prepare AD4Dpf', 'Macro')

adtlib.addNode(Assign_AD4Types, 'AD4_typer', 'Mapper')
adtlib.addNode(Add_SolvationParameters, 'Add Solvation Parameters', 'Mapper')
adtlib.addNode(ManageRotatableBonds, 'Manage Rotatable Bonds', 'Mapper')
adtlib.addNode(MergeNonPolarHydrogens, 'Merge NonPolar Hydrogens', 'Mapper')
adtlib.addNode(MergeLonePairs, 'Merge Lone Pairs', 'Mapper')
adtlib.addNode(RemoveWaters, 'Remove Water Residues', 'Mapper')
adtlib.addNode(ManageAromaticCarbons, 'Manage Aromatic Carbons', 'Mapper')

adtlib.addNode(Ligand_Writer, 'Ligand Writer', 'Output')
adtlib.addNode(Ligand_Writer_AD4, 'AD4 Ligand Writer', 'Output')


try:
    UserLibBuild.addTypes(adtlib, 'MolKit.VisionInterface.MolKitTypes')
except Exception, e:
    print "loading types failed:", e

