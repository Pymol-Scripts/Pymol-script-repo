## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#########################################################################
#
# Date: Nov 2001 Authors: Michel Sanner, Daniel Stoffler
#
#    sanner@scripps.edu
#    stoffler@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner, Daniel Stoffler and TSRI
#
#########################################################################

import warnings
import numpy.oldnumeric as Numeric, string, os

from NetworkEditor.items import NetworkNode
from Vision import UserLibBuild

from DejaVu.VisionInterface.GeometryNodes import GeometryNode

from MolKit.molecule import Atom, AtomSet, Molecule, MoleculeSet
from MolKit.protein import Residue, ResidueSet, Chain, ChainSet
from MolKit.moleculeWriter import MoleculeWriter
from MolKit.pdbWriter import PdbWriter, PdbqWriter, PdbqsWriter, PdbqtWriter
from MolKit.pqrWriter import PqrWriter
from MolKit.pdbParser import PdbParser
from MolKit.chargeCalculator import GasteigerChargeCalculator, KollmanChargeCalculator
from MolKit.hydrogenBuilder import HydrogenBuilder, PolarHydrogenBuilder

try:
    from mslib import MSMS, _msms
    msmsFound=1
except ImportError:
    import traceback
    traceback.print_exc()
    msmsFound=0


def importMolKitLib(net):
    try:
        from MolKit.VisionInterface.MolKitNodes import molkitlib
        net.editor.addLibraryInstance(
            molkitlib, 'MolKit.VisionInterface.MolKitNodes', 'molkitlib')
    except:
        import traceback
        traceback.print_exc()
        warnings.warn(
            'Warning! Could not import molitlib from MolKit.VisionInterface')


def importVizLib(net):
    try:
        from DejaVu.VisionInterface.DejaVuNodes import vizlib
        net.editor.addLibraryInstance(
            vizlib, 'DejaVu.VisionInterface.DejaVuNodes', 'vizlib')
    except:
        import traceback
        traceback.print_exc()
        warnings.warn(
            'Warning! Could not import vizlib from DejaVu/VisionInterface')


def importSymServLib(net):
    try:
        from symserv.VisionInterface.SymservNodes import symlib
        net.editor.addLibraryInstance(
            symlib, 'symserv.VisionInterface.SymservNodes', 'symlib')
    except:
        import traceback
        traceback.print_exc()
        warnings.warn(
            'Warning! Could not import symlib from symserv/VisionInterface')


class writePDBFiles(NetworkNode):
    """
    Write a PDF file for each set of coordinates in coordsList

    [filenames] <- writePDBFiles(mol, coordsList, filename)

    - mol is a MolKit molecule object
    - coordsList is a list of atomic coordinates sets. Each set has to match
      the number of atoms in mol or be a multiple of the number of atoms in mol.
      In the latter case the list will be split.
    - filename is a string used to generate name o fthe files by concatenating
      of zero padded integer ans .pdb
    """

    mRequiredTypes = {'MoleculeType': 'MolKit.VisionInterface.MolKitTypes'}
    mRequiredSynonyms = [
    ]
    def __init__(self, constrkw = {},  name='writePDBFiles', **kw):
        kw['constrkw'] = constrkw
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)
        code = """from MolKit.pdbWriter import PdbWriter
import numpy

def doit(self, mol, coordsList, filename):
    li = len(coordsList[0])
    allAtoms = mol.allAtoms
    la = len(allAtoms)
    if li != la:
        li = len(coordsList)
        assert li % la == 0
        coordsList = numpy.array(coordsList)
        coordsList.shape = (li/la, la, 3)

    origCoords = allAtoms.coords
    fnames = []
    try:
        for i, coords in enumerate(coordsList):
            fname = filename+'''%05d'''%i
            fnames.append(fname)
            allAtoms.coords = coords
            writer = PdbWriter()
            writer.write(filename+'''%05d'''%i, mol)
    finally:
        allAtoms.coords = origCoords
	pass

    self.outputData(filenames=fnames)
"""
        self.configure(function=code)

        ip = self.inputPortsDescr
        ip.append( {'name': 'mol', 'datatype': 'Molecule'} )
        ip.append( {'name': 'coordsList', 'datatype': 'list'} )
        ip.append( {'name': 'filename', 'datatype': 'string'} )

        self.outputPortsDescr.append(
            {'name': 'filenames', 'datatype': 'list'} )

        self.widgetDescr['filename'] = {
            'initialValue': 'test',
            'labelGridCfg': {'column': 0, 'row': 0},
            'master': 'node',
            'widgetGridCfg': {'labelSide': 'left', 'column': 1, 'row': 0},
            'labelCfg': {'text': 'filename:'}, 'class': 'NEEntry'}


    def beforeAddingToNetwork(self, net):
        try:
            ed = net.getEditor()
        except:
            import traceback; traceback.print_exc()
            print 'Warning! Could not import widgets'


class PDB_SYMTRY(NetworkNode):
    """inputs PDB file, parses SYMTRY records and returns a list of (4x4)
    matrices."""

    def __init__(self, name='PDB_SYMTRY', **kw):
        kw['name']=name
        apply( NetworkNode.__init__, (self,), kw)

	code = """def doit(self, mol):
    if mol:
        matrices = []
        all = lines = mol.data[0].parser.allLines
        rem290 = filter(lambda x: x[:18]=="REMARK 290   SMTRY", all)
        for i in range(0, len(rem290), 3):
            lines = []
            mat = Numeric.identity(4).astype('f')
            for j in range(3):
               line1 = rem290[i+j]
               r1 = float(line1[24:33])
               r2 = float(line1[34:43])
               r3 = float(line1[44:53])
               t = float(line1[54:68])
               mat[j] = (r1, r2, r3, t)
            mat[3] = (0., 0., 0., 1.)
            #matrices.append( Numeric.array( Numeric.transpose( mat) ) )
            matrices.append( mat )

        self.outputData(matrices=matrices)\n"""

        if code: self.setFunction(code)

	ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecule')
        
	op = self.outputPortsDescr
        op.append(datatype='instancemat(0)', name='matrices')


    def beforeAddingToNetwork(self, net):
        # import molkitlib
        importSymServLib(net)
           


class PDB_MTRIX(NetworkNode):
    """
    Get non crystallographic symmetry matrices

    inputs PDB file, parses MTRIX1,2,3 records and output a list of (4x4)
    instance matrices.
    """

    def __init__(self, name='PDB_MTRIX', **kw):
        kw['name']=name
        apply( NetworkNode.__init__, (self,), kw)

	code = """def doit(self, mol):
    if mol:
        all = lines = mol.data[0].parser.allLines
        mtrix = filter(lambda x: x[:5]=="MTRIX", all)
        import numpy
        matrices = numpy.zeros( (len(mtrix)/3, 4, 4), 'f')
        for i in range(0, len(mtrix), 3):
            lines = []
            mat = Numeric.identity(4).astype('f')
            for j in range(3):
               line1 = mtrix[i+j]
               num = int(line1[7:10])-1
               rx = float(line1[10:20])
               ry = float(line1[20:30])
               rz = float(line1[30:40])
               t = float(line1[45:55])
               matrices[num,j] = (rx, ry, rz, t)
            matrices[num,3,3] = 1.

        mat = []
        for m in matrices:
            mat.append(m)
        self.outputData(matrices=mat)\n"""

        if code: self.setFunction(code)

	ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecule')
        
	op = self.outputPortsDescr
        op.append(datatype='instancemat(0)', name='matrices')


    def beforeAddingToNetwork(self, net):
        # import molkitlib
        importSymServLib(net)
           


class CenterAndRadius(NetworkNode):

    def __init__(self, name='centerAndRad', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')

        op = self.outputPortsDescr
        op.append(datatype='None', name='center')
        op.append(datatype='float', name='radius')

        code = """def doit(self, molecules):
    import numpy.oldnumeric as Numeric
    from math import sqrt
    centers = []
    radii = []
    for mol in molecules:
        c = Numeric.array(mol.chains.residues.atoms.coords)
        g = Numeric.add.reduce(c) / len(c)
        centers.append( g )
        d = c-[g]
        radii.append( sqrt(max(Numeric.add.reduce(d*d, 1))) )

    if len(centers):
        self.outputData(center=centers, radius=radii)
"""
        if code: self.setFunction(code)



class ReadMolecule(NetworkNode):
    """based on the MolKit.Read fucntion. Reads any file format recognized
by MolKit. This currently includes PDB, PDBQ, PDBQS, MOL2 and PQR.
Parse-able keywords can be specified in the param.panel of this node.
Input:  filename (string)
Output: MoleculeSet"""
    
    def __init__(self, name='Read Molecule', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        fileTypes=[('pdb', '*.pdb'),('pdbq', '*.pdbq'),('pdbqs', '*.pdbqs'),
                   ('pdbqt', '*.pdbqt'),('pqr', '*.pqr'),('mol2', '*.mol2'),('all', '*')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'read molecule', 'width':10,
            'labelCfg':{'text':'file:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='filename')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='MolSets')

        code = """def doit(self, filename):
    if filename and len(filename):
        from MolKit import Read
        mols = Read(filename)
        if mols:
            self.outputData(MolSets=mols)\n"""

        self.setFunction(code)

    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass
    

class WritePDB(NetworkNode):
    """Save as PDB."""
    
    def __init__(self, name='Write PDB', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        fileTypes=[('pdb', '*.pdb'), ('all', '*')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileSaver', 'master':'node',
            'filetypes':fileTypes, 'title':'save PDB', 'width':10,
            'labelCfg':{'text':'file:'},
            }

        self.inputPortsDescr.append(datatype='list', name='nodes')
        self.inputPortsDescr.append(datatype='string', name='filename')

        code = """def doit(self, nodes, filename):
    if filename and len(filename):
        mol = nodes.top.uniq()
        assert len(mol)==1

        recType = 'all'  # this could be exposed to the user?
        sort = 1

        if sort:
            nodes.sort()

        mol[0].write(filename, PdbWriter(), nodes, recType )\n"""

        self.setFunction(code)


class WriteMolecule(NetworkNode):
    """Write any file format recognized
by MolKit. This currently includes PDB, PDBQ, PDBQS, MOL2 and PQR.
Parse-able keywords can be specified in the param.panel of this node.
Input: MoleculeSet
Input:  filename (string)"""
    
    def __init__(self, name='Write Molecule', **kw):
        self.writers = {}
        self.writers['pdb'] = PdbWriter
        self.writers['pqr'] = PqrWriter
        self.writers['pdbq'] = PdbqWriter
        self.writers['pdbqs'] = PdbqsWriter
        self.writers['pdbqt'] = PdbqtWriter
        self.required_attrs = {}
        self.required_attrs['pdb'] = [] #?
        self.required_attrs['pqr'] = ['charge', 'radius'] #same as pqrRadius
        self.required_attrs['pdbq'] = ['charge']
        self.required_attrs['pdbqs'] = ['charge', 'AtSolPar', 'AtVol']
        self.required_attrs['pdbqt'] = ['charge', 'autodock_element']

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        fileTypes=[('pdb', '*.pdb'),('pdbq', '*.pdbq'),('pdbqs', '*.pdbqs'),
                   ('pdbqt', '*.pdbqt'),('pqr', '*.pqr')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'write molecule', 'width':10,
            'labelCfg':{'text':'output file:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')
        ip.append(datatype='string',  name='filename')
        ip.append(datatype='int', required=False, name='sort', defaultValue=False)
        ip.append(datatype='list', required=False, name='records', defaultValue=['ATOM','HETATM'])

        op = self.outputPortsDescr
        op.append(datatype='string', name='filename')

        code = """def doit(self, mols, filename, sort, records):
    if mols and filename and len(filename):
        ext = filename.split('.')[-1]
        if not self.writers.has_key(ext):
            print 'unrecognized outputfile type', ext
            return 
        writer = self.writers[ext]()
        #check for required attributes
        #pqr requires charges, radius
        #pdbq requires charges
        #pdbqt requires charges and types
        #pdbqs requires charges and solvation parameters
        mol = mols[0]
        ats = mol.allAtoms
        for attr in self.required_attrs[ext]:
            try:
                getattr(ats, attr)
            except:
                print "all atoms do not have required attribute:", attr
                return 'ERROR'
        #need support for specifying records (?) 
        #NB: pqrWriter doesnot have same signature: no records parameter
        try:                     #should work if ext!='pqr'
            writer.write(filename, mol.allAtoms, sort=sort, records=records)
        except:
            writer.write(filename, mol.allAtoms, sort=sort)
        self.outputData(filename=filename)\n"""

        self.setFunction(code)

    def myCallback(self, event=None):
        #self.paramPanel.run()
        pass
    

class AssignRadii(NetworkNode):

    def __init__(self, name='Assign Radii', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['united'] = {
            'class':'NECheckButton', 'master':'node',
            'labelCfg':{'text':'united radii'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')
        ip.append(datatype='int', name='united')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='molecules')
        op.append(datatype='list', name='radii')
        
        code = """def doit(self, molecules, united ):
    rad = []
    for mol in molecules:
        rad.extend(mol.defaultRadii(united=united))
    if rad:
        self.outputData(molecules=molecules, radii=rad)
"""

        self.setFunction(code)


class CalculateGasteigerCharges(NetworkNode):

    def __init__(self, name='Calculate Gasteiger Charges', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='mols')
        op.append(datatype='float', name='charge_total')
        
        code = """def doit(self, mols ):
    calculator = GasteigerChargeCalculator()
    for mol in mols:
        calculator.addCharges(mol.allAtoms)
    import numpy.oldnumeric as Numeric
    charge_total = Numeric.add.reduce(mols.allAtoms.charge)
    self.outputData(mols=mols, charge_total=charge_total)
"""

        self.setFunction(code)


class AddKollmanCharges(NetworkNode):

    def __init__(self, name='Add Kollman Charges', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='mols')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='mols')
        op.append(datatype='float', name='charge_total')
        
        code = """def doit(self, mols ):
    calculator = KollmanChargeCalculator()
    for mol in mols:
        calculator.addCharges(mol.allAtoms)
    import numpy.oldnumeric as Numeric
    charge_total = Numeric.add.reduce(mols.allAtoms.charge)
    self.outputData(mols=mols, charge_total=charge_total)
"""

        self.setFunction(code)


class AddHydrogens(NetworkNode):

    def __init__(self, name='Add Hydrogens', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='molecules', balloon='molecules with added hydrogens')
        op.append(datatype='int', name='new_h_ct', balloon='number of hydrogens added to molecules')
        
        code = """def doit(self, molecules):
    h_builder = HydrogenBuilder()
    for mol in molecules:
        new_h_ct = h_builder.addHydrogens(mol)
    self.outputData(molecules=molecules, new_h_ct=new_h_ct)
"""

        self.setFunction(code)


class AddPolarHydrogens(NetworkNode):

    def __init__(self, name='Add Polar Hydrogens', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='molecules')
        op.append(datatype='int', name='new_h_ct')
        
        code = """def doit(self, molecules):
    h_builder = PolarHydrogenBuilder()
    for mol in molecules:
        new_h_ct = h_builder.addHydrogens(mol)
    self.outputData(molecules=molecules, new_h_ct=new_h_ct )
"""

        self.setFunction(code)


class ParseCryst1(NetworkNode):

    def __init__(self, name='ParseCryst1', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='molecules')

        code = """def doit(self, molecules):
    for mol in molecules:
        p = mol.parser
        if not isinstance(p, PdbParser):
            print 'WARNING molecule %s was not read from a PDB file'%mol.name
            continue
        p.parse_PDB_CRYST1( p.getRecords(p.allLines, 'CRYST1') )

    self.outputData(molecules=molecules)
"""

        self.setFunction(code)
        

class NodeSelector(NetworkNode):
    """This node selects subsets of nodes (i.e. Atoms, Residues, Chains, ...)
    from a list of incomming molecules.
    First all nodes of the currently selected level are found, then a
    selection operation is performed on these nodes.
    Selection can be specified using regular expressions or lambda functions.

    Example:
          # set the level to Residue and NOTE that the port name changes :)
          type ^G* to select all residues with names starting with a 'G'
          type lambda x: x._uniqIndex < 10 to select the 10 first residues
"""
    def setLevel(self, level):
        if len(self.outputPorts)==0:
            return # we got here before node is fully created
        newlevel = self.levels[level][0]
        if newlevel==self.level:
            return
        self.level = newlevel
        p = self.outputPorts[0]
        p.setDataType(self.levels[level][1].__name__)

        
    def __init__(self, name='Select Nodes', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        #self.readOnly = 1
        self.levels = {'Atom': (Atom, AtomSet),
                       'Residue': (Residue, ResidueSet),
                       'Chain': (Chain, ChainSet),
                       'Molecule': (Molecule, MoleculeSet) }

        self.level = Atom

        self.widgetDescr['selectionString'] = {
            'class':'NEEntry', 'master':'node', 'width':14, 
            'labelGridCfg':{'sticky':'w'},
            'widgetGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'select:'},
            }

        self.widgetDescr['nodeType'] = {
            'class':'NEComboBox', 'master':'node',
            'choices':self.levels.keys(),
            'fixedChoices':True,
            'initialValue':'Atom',
            'entryfield_entry_width':8,
            'selectioncommand':self.setLevel,
            'labelGridCfg':{'sticky':'w'},
            'widgetGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'level:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='TreeNodeSet', name='nodes')
        ip.append(datatype='string', name='nodeType')
        ip.append(datatype='string', required=False, name='selectionString')

        op = self.outputPortsDescr
        op.append(datatype='AtomSet', name='nodes')
        
        code = """def doit(self, molecules,  nodeType, selectionString):
    elemType, setType = self.levels[nodeType]
    result = setType([])
    for mol in molecules:
        selnodes = setType(mol.findType(elemType))
        if selectionString:
            selnodes = selnodes.get(selectionString)
        result = result + selnodes
    if len(result):
        result.sort()
        self.outputData(nodes=result)\n"""

        self.setFunction(code)

class RMSDFromTwoAtomSet(NetworkNode):

    def __init__(self, name='RMSD', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atomSet1')
        ip.append(datatype='AtomSet', name='atomSet2')

        op = self.outputPortsDescr
        op.append(datatype='float', name='rmsd')

        from mglutil.math.rmsd import RMSDCalculator
        self.RMSD = RMSDCalculator()

        code = """def doit(self, atomSet1, atomSet2):
    if atomSet1 is None or atomSet2 is None:
        return
    if len(atomSet1)!=len(atomSet2):
        return 

    # sort the two atomSet
    needSorting = False
    for i in range(len(atomSet1)):
        if atomSet1[i].name != atomSet2[i].name:
            needSorting = True
            break

    
    if needSorting:
        tmp = AtomSet()
        names = atomSet2.name
        for atom in atomSet1:
            try:
                idx = names.index(atom.name)
                tmp.append(atomSet2[idx])
            except:
                print atom.name, 'not found in', atomSet2
        mob = tmp.coords
    else:
        mob = atomSet2.coords
        
    ref = atomSet1.coords
    self.RMSD.setRefCoords(ref)
    rmsd = self.RMSD.computeRMSD(mob)
    
    self.outputData(rmsd=rmsd)

"""

        self.setFunction(code)


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)
        

class AtomsAsCPK(GeometryNode):

    def __init__(self, name='CPK', **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'CPK'), kw )

        #self.readOnly = 1

#        from DejaVu.Spheres import Spheres
#        self.cpk = Spheres('CPK', inheritMaterial=0)

        self.widgetDescr['quality'] = {
            'class':'NEThumbWheel', 'master':'node',
            'width':60, 'height':18, 'oneTurn':10,
            'lockType':1, 'min':3, 'increment':1, 'type':'int', 'wheelPad':2,
            'initialValue':8,
            'labelCfg':{'text':'quality:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atoms')
        
        ip.append(datatype='float(0)', required=False, name='radii')
        #ip.append(datatype='vector', required=False, name='radii')
        
        ip.append(datatype='colorfloat3or4(0)', required=False, name='colors')
        ip.append(datatype='int', required=False, name='quality')


        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        code = """def doit(self, atoms, radii, colors, quality,
name, geoms, instanceMatrices, geomOptions, parent):
    GeometryNode.doit(self, name, geoms, instanceMatrices, geomOptions, parent)
    if atoms is None:
        return
    try:
        len(radii)
        islist = True
    except TypeError:
        islist = False
    if islist:
        assert len(radii)==len(atoms.coords)
    if colors is None:
        colors = ( (0.7, 0.7, 0.7), )
        
        
    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=atoms.coords, radii=radii, materials=colors, 
              quality=quality, inheritMaterial=inheritMaterial)
        self.outputData(CPK=g, allGeometries=self.geoms)
    else:
        self.outputData(CPK=None, allGeometries=self.geoms)

"""
        self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu.Spheres import Spheres
        self.geoms.append(Spheres(name))
        return 1 # num of appended geoms


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)
        

class AtomsAsSticks(NetworkNode):

    def __init__(self, name='Sticks', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.indices     = None

        #self.readOnly = 1

        from DejaVu.Cylinders import Cylinders
        self.sticks = Cylinders('sticks', inheritMaterial=0)

        self.widgetDescr['cradii'] = {
            'class':'NEThumbWheel', 'master':'node',
            'width':60, 'height':18, 'wheelPad':2, 'oneTurn':1,
            'lockType':1, 'min':0, 'increment':.1, 'type':'float',
            'initialValue':0.25,
            'labelCfg':{'text':'radius:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atoms')
        ip.append(datatype='colorsRGB', required=False, name='colors')
        ip.append(datatype='instancemat(0)', required=False,
                  name='instance matrices')
        ip.append(datatype='float', required=False, name='cradii')

        op = self.outputPortsDescr
        op.append(datatype='geom', name='sticks')

        code = """def doit(self, atoms, colors, instanceMatrices, cradii):
    if atoms is None:
        return
    if colors is None:
        colors = ( (0.7, 0.7, 0.7), )

    if instanceMatrices:
        mat = Numeric.array(instanceMatrices)
    else:
        mat = [Numeric.identity(4).astype('f')]

    if cradii is None:
        cradii = 0.25
    bonds, atnobnd  = atoms.bonds
    self.indices = map(lambda x: (x.atom1._bndIndex_,
                                  x.atom2._bndIndex_), bonds)
    self.sticks.Set(vertices=atoms.coords, faces=self.indices,
                    radii=cradii, materials=colors)
    self.sticks.Set(instanceMatrices=mat)

    self.outputData(sticks=self.sticks)\n"""

        self.setFunction(code)


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)


class SESVertices(NetworkNode):
    """Extract specific vertices from the SES surface

Input ports:
    MSMS: an MSMSobject instance as output by the MSMS calculation node
    nodes: a TreeNodeSet describing parts of a molecule. This will be
           expanded to an AtomSet triangles for which at least 'selnum'
           vertices belong to atoms in this set will be extracted.
           If ommited, the whole surface is extracted.
    selnum: number of vertices per triangle that need to belong to atoms
            in nodes for a triangle to be output.
    component: MSMS component (i.e. a closed surface) to be extracted.
               Component 0 is always the external component.
    ftype: can be any of 'contact', 'reentrant' or 'toric' to select only
           vertices in this type of faces
    vtype: can be any of 'inside', 'edge' or 'vertices' to select only
           vertices inside contact face, on edges or vertices respectively
           
OutputPort:
    vertices: array of X,Y,Z coordinates
    vnormals: array of nx,ny,nz vertices normal vectors
"""
    def __init__(self, name='gSESVertices', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['selnum'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':3, 'min':1, 'max':3, 'lockBMax': 1,
            'lockMin': 1, 'type':'int',
            'initialValue':3,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'nbVert'},
            }
        
        self.widgetDescr['component'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':5, 'min':0, 'lockMin':1, 'lockMax':1, 'type':'int',
            'initialValue':0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'component'},
            }

        self.widgetDescr['ftype'] = { 'class': 'NEComboBox',
            'initialValue': 'contact', 'fixedChoices': 1,
            'entryfield_entry_width':12,
            'choices': ['contact', 'reentrant', 'toric'],
            'labelGridCfg':{}, 'master': 'node', 'widgetGridCfg': {},
            'labelCfg': {'text': 'vtype:'},
            }

        self.widgetDescr['vtype'] = { 'class': 'NEComboBox',
            'initialValue': 'inside', 'fixedChoices': 1,
            'entryfield_entry_width':12,
            'choices': ['inside', 'edges', 'vertices'],
            'labelGridCfg':{}, 'master': 'node', 'widgetGridCfg': {},
            'labelCfg': {'text': 'vtype:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='MSMSobject', name='MSMS')
        ip.append(datatype='TreeNodeSet',
                  balloon='Atoms for which triangles are extracted',
                  required=False, name='nodes')
        ip.append(datatype='int', required=False, name='selnum', defaultValue=3)
        ip.append(datatype='int', required=False, name='component', defaultValue=0)
        ip.append(datatype='string', required=False, name='ftype', defaultValue='contact')
        ip.append(datatype='string', required=False, name='vtype', defaultValue='inside')

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='SAScenters')
        op.append(datatype='normals3D', name='normals')

        code = """def doit(self, msms, nodes, selnum, component, ftype, vtype):
    vf,vi,f = msms.getTriangles()
    v = {}
    n = {}
    if ftype=='contact': faceTest=1
    elif ftype=='reentrant': faceTest=2
    elif ftype=='toric': faceTest=3
    else: raise(ValueError, 'bad face type')
    
    for fa in f:
        if fa[3]==faceTest:
            for i in range(3):
                v1 = fa[i]
                vt = vi[v1][0]
                 # pick vertices based on vtype
                if (vtype=='inside' and vt>0) or \
                   (vtype=='edges' and vt<0) or \
                   (vtype=='vertices' and vt<0):
                    v[v1] = vf[v1][0:3]
                    n[v1] = vf[v1][3:6]
    self.outputData(SAScenters=v.values())
    self.outputData(normals=n.values())
"""

        self.setFunction(code)


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)

from Pmv.msmsParser import MSMSParser
class ReadMSMS(NetworkNode):
    """based on the MSMSParser object. Reads an MSMS surface
Input:
    filename (string), filename.vert and filename.face will be read
Output:
    vertices: array of vertices
    normals:  array of normals
    vprop:    array of vertex properties (integer values)
    tri:      array of faces
    triprop:  array of face properties (integers)"""
    
    def __init__(self, name='Read MSMS', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        fileTypes=[('msms', '*.vert')]

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node',
            'filetypes':fileTypes, 'title':'read molecule', 'width':10,
            'labelCfg':{'text':'file:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='filename')

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='vertices')
        op.append(datatype='faceIndices', name='tri')
        op.append(datatype='normals3D', name='normals')
        op.append(datatype='int(0,3)', name='vprop')
        op.append(datatype='int(0,2)', name='triprop')

        code = """def doit(self, filename):
    if filename and len(filename):
        msmsParser = MSMSParser()
        msmsParser.parse(filename, filename[:-4]+'face')
        self.outputData(
            vertices=msmsParser.vertices,
            tri=msmsParser.faces,
            normals=msmsParser.normals,
            vprop=msmsParser.vertProp ,
            triprop=msmsParser.facesProp)"""

        self.setFunction(code)


try:
    bhtreelibFound = True
    from bhtree import bhtreelib

    class SurfaceAtoms(NetworkNode):
        """outputs a list of atoms contributing to the MSMS serface

    Input ports:
        MSMS: an MSMSobject instance as output by the MSMS calculation node
        nodes: a TreeNodeSet describing parts of a molecule. This will be
               expanded to an AtomSet triangles for which at least 'selnum'
               vertices belong to atoms in this set will be extracted.
               If ommited, the whole surface is extracted.
        component: MSMS component (i.e. a closed surface) to be extracted.
                   Component 0 is always the external component.
        cut0ff: if an atom is less than cutoff away from the surface it will be
                selected

    OutputPort:
        surfaceAtoms: AtomSet of atoms close to the surface
        interiorAtoms: AtomSet of atoms not close to the surface
    """
        def __init__(self, name='SurfaceAtoms', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw )

            self.widgetDescr['component'] = {
                'class':'NEDial', 'size':50,
                'oneTurn':5, 'min':0, 'lockMin':1, 'lockMax':1, 'type':'int',
                'initialValue':0,
                'labelGridCfg':{'sticky':'w'},
                'labelCfg':{'text':'component'},
                }
            self.widgetDescr['cutoff'] = {
                'class':'NEDial', 'size':50,
                'oneTurn':5.0, 'min':0, 'lockMin':1, 'type':'float',
                'initialValue':0.0,
                'labelGridCfg':{'sticky':'w'},
                'labelCfg':{'text':'cutoff'},
                }

            ip = self.inputPortsDescr
            ip.append(datatype='MSMSobject', name='MSMS')
            ip.append(datatype='TreeNodeSet',
                      balloon='Atoms for which triangles are extracted',
                      required=False, name='nodes')
            ip.append(datatype='int', required=False, name='component', defaultValue=0)
            ip.append(datatype='float', required=False, name='cutoff', defaultValue=0.)

            op = self.outputPortsDescr
            op.append(datatype='AtomSet', name='surfaceAtoms')
            op.append(datatype='AtomSet', name='interiorAtoms')


            code = """def doit(self, msms, nodes, component, cutoff):
        vf,vi,f = msms.getTriangles()
        from mglutil.util.uniq import uniq
        atInd = uniq(vi[:,1])
        atoms = msms.atoms
        surfats ={}
        intats = []

        if cutoff is None: cutoff=0.0
        
        # find surface atoms
        for n in atInd:
            surfats[atoms[n-1]]=n-1

        # find complementary set
        for a in atoms:
            if surfats.get(a, None) is None:
                intats.append(a)

        if cutoff!=0.0:
            import numpy.oldnumeric as Numeric
            points = Numeric.array(vf[:,:3], 'f')
            ids = Numeric.arange(len(points)).astype('i')
            result = Numeric.zeros( (len(points),) ).astype('i')
            bht = bhtreelib.TBHTree( points, ids, 10, 10, 9999.0 )
            for a in intats:
                nb = bht.ClosePoints( tuple(a.coords), cutoff, result)
                if nb:
                    print 'atom' , a, 'is close'
                    surfats[a] = True
                    intats.remove(a)

        self.outputData(surfaceAtoms=AtomSet(surfats.keys()))
        self.outputData(interiorAtoms=AtomSet(intats))
"""

            self.setFunction(code)

except:
      bhtreelibFound = False


class AtomsAsMSMS(NetworkNode):
    """Compute solvent excluded surface using the MSMS library.

    Input ports:
      atoms: a set of atoms for which to compute a surface
      radii: optional atomic radii. If they are ommited, we tri to get them
             from the atoms
      seedAtoms: a list of up to 3 atoms that are used to position the rolling
             probe initially. Can be used to compute cavity surfaces.
      density: floating point value describing the number of vertices per
             angstrom square used for triangulating the surface
      probe_radius: floating point value decribing the radius of the rolling
             probe sphere.
      allComp: Boolean value triggering the computation of all components
      
      OutputPort:
        MSMS: outputs an MSMSObject

      The set of atoms used to calculate the surface are stored in the atoms
      attribute of the MSMSobject. By default, the probe_radius and the
      density parameters are bound to dials and allComp is bound to a check-
      button widget located in the parameter panel.
      """
    def __init__(self, name='MSMS', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        #self.readOnly = 1

        from DejaVu.IndexedPolygons import IndexedPolygons
        self.srf = IndexedPolygons('MSMS', inheritMaterial=0)

        self.widgetDescr['density'] = {
            'class':'NEDial', 'size':50, 
            'oneTurn':10, 'min':1.0,
            'increment':0.1, 'type':'float',
            'initialValue':1.0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'density'},
            }
        
        self.widgetDescr['probe_radius'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':10, 'min':0.5, 'increment':0.1, 'type':'float',
            'initialValue':1.5,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'probe radius'},
            }
        
        self.widgetDescr['allComp'] = {
            'class':'NECheckButton',
            'initialValue':0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'all Components'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atoms')
        ip.append(datatype='vector', required=False, name='radii')
        ip.append(datatype='AtomSet', required=False, name='seedAtoms')
        ip.append(datatype='float', required=False, name='density', defaultValue=1.)
        ip.append(datatype='float', required=False, name='probe_radius', defaultValue=1.5)
        ip.append(datatype='int', required=False, name='allComp', defaultValue=0)
        
        op = self.outputPortsDescr
        op.append(datatype='MSMSobject', name='MSMS')
        
        code = """def doit(self, atoms, radii, seedAtoms, density, probe_radius, allComp):
    if atoms is None:
        return
    if allComp is None:
        allComp = 0
    if radii is None:
        try:
            radii = atoms.radius
        except:
            print 'No Radii'
            return

    if max(radii)==0.0:
        return 'STOP'
        
    # create MSMS object
    srf = self.msms = MSMS(coords=atoms.coords, radii=radii )
    srf.atoms = atoms
    srf.all_components = allComp

    # if atoms to be used to seed the probe a specified, use them
    if seedAtoms is not None:
        seedInd = [-1,-1,-1]
        i=0
        for a in seedAtoms:
            seedInd[i] = atoms.index(a)
            i = i+1
        srf.rsr.ffnba = len(seedAtoms)
        srf.rsr.set_ffat(seedInd[0], seedInd[1], seedInd[2])

    # calculate the surface
    srf.compute(probe_radius=probe_radius, density=density)

    self.outputData(MSMS=srf)
"""

        self.setFunction(code)



class SpheresAsMSMS(NetworkNode):
    """Compute solvent excluded surface using the MSMS library.

    Input ports:
      centers: a set of atoms for which to compute a surface
      radii: optional atomic radii. If they are ommited, we tri to get them
             from the atoms
      seedAtoms: a list of up to 3 atoms that are used to position the rolling
             probe initially. Can be used to compute cavity surfaces.
      density: floating point value describing the number of vertices per
             angstrom square used for triangulating the surface
      probe_radius: floating point value decribing the radius of the rolling
             probe sphere.
      allComp: Boolean value triggering the computation of all components
      
      OutputPort:
        MSMS: outputs an MSMSObject

      The set of atoms used to calculate the surface are stored in the atoms
      attribute of the MSMSobject. By default, the probe_radius and the
      density parameters are bound to dials and allComp is bound to a check-
      button widget located in the parameter panel.
      """
    def __init__(self, name='MSMS from spheres', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        #self.readOnly = 1

        from DejaVu.IndexedPolygons import IndexedPolygons
        self.srf = IndexedPolygons('MSMS', inheritMaterial=0)

        self.widgetDescr['density'] = {
            'class':'NEDial', 'size':50, 
            'oneTurn':10, 'min':1.0,
            'increment':0.1, 'type':'float',
            'initialValue':1.0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'density'},
            }
        
        self.widgetDescr['probe_radius'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':10, 'min':0.5, 'increment':0.1, 'type':'float',
            'initialValue':1.5,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'probe radius'},
            }
        
        self.widgetDescr['allComp'] = {
            'class':'NECheckButton',
            'initialValue':0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'all Components'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='coordinates3D', name='centers')
        ip.append(datatype='None', name='radii')
        ip.append(datatype='list', required=False, name='seedIndices')
        ip.append(datatype='float', name='density')
        ip.append(datatype='float', name='probe_radius')
        ip.append(datatype='int', name='allComp')
        
        op = self.outputPortsDescr
        op.append(datatype='MSMSobject', name='MSMS')
        
        code = """def doit(self, centers, radii, seedIndices, density, probe_radius, allComp):

    if allComp == None:
        allComp = 0

    try:
        len(radii)
    except:
        radii= (radii,)*len(centers)

    # create MSMS object
    srf = self.msms = MSMS(coords=centers, radii=radii )
    srf.all_components = allComp

    # if atoms to be used to seed the probe a specified, use them
    if seedIndices is not None:
        seedInd = [-1,-1,-1]
        seedInd[:len(seedIndices)] = seedIndices
        srf.rsr.ffnba = len(seedIndices)
        srf.rsr.set_ffat(seedInd[0], seedInd[1], seedInd[2])

    # calculate the surface
    srf.compute(probe_radius=probe_radius, density=density)

    self.outputData(MSMS=srf)
"""

        self.setFunction(code)


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)


class SaveMSMS(NetworkNode):
    """save the triangulated Solvent Excluded Surface in ascii format"""
    
    def __init__(self, constrkw = {},  name='save MSMS', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr
        ip.append(datatype='MSMSobject', name='MSMS')
        ip.append(datatype='string', name='filename')

        code = """def doit(self, MSMS, filename):
    MSMS.write_triangulation(filename)
"""
        if code: self.setFunction(code)


class GetMSMSareas(NetworkNode):
    """Compute and output surface areas MSMS component.

Input ports:
    MSMS: an MSMSobject instance as output by the MSMS calculation node
    component: MSMS component (i.e. a closed surface) to be extracted.
               Component 0 is always the external component.
OutputPort:
    ses_area: surface area of the solvent excluded surface
    sas_area: surface area of the solvent accessible surface
    """
    def __init__(self, name='get MSMS areas', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['component'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':5, 'min':0, 'lockMin':1, 'lockMax':1, 'type':'int',
            'initialValue':0,
            'labelCfg':{'text':'component'},
            }
        
        ip = self.inputPortsDescr
        ip.append(datatype='MSMSobject', name='MSMS')
        ip.append(datatype='int', required=False, name='component', defaultValue=0)

        op = self.outputPortsDescr
        op.append(datatype='float', name='ses_area')
        op.append(datatype='float', name='sas_area')

        code = """def doit(self, msms, component):
    self.inputPorts[1].widget.configure(max=msms.rsr.nb-1, type='int')
    msms.compute_ses_area()
    ses = msms._getSESComp(component)
    self.outputData(ses_area=ses.a_ses_area)
    self.outputData(sas_area=ses.a_sas_area)\n"""

        self.setFunction(code)
        

class GetMSMStriangles(NetworkNode):
    """Extract triangles from a given MSMS component and for a given set of
atoms.

Input ports:
    MSMS: an MSMSobject instance as output by the MSMS calculation node
    nodes: a TreeNodeSet describing parts of a molecule. This will be
           expanded to an AtomSet triangles for which at least 'selnum'
           vertices belong to atoms in this set will be extracted.
           If ommited, the whole surface is extracted.
    selnum: number of vertices per triangle that need to belong to atoms
            in nodes for a triangle to be output.
    component: MSMS component (i.e. a closed surface) to be extracted.
               Component 0 is always the external component.
    
OutputPort:
    vertices: array of X,Y,Z coordinates
    indices:  array of 0-based i,j,k indices describing triangular facets
    vnormals: array of nx,ny,nz vertices normal vectors
"""
    def __init__(self, name='get MSMS triangles', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['selnum'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':3, 'min':1, 'max':3, 'lockBMax': 1,
            'lockMin': 1, 'type':'int',
            'initialValue':3,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'nbVert'},
            }
        
        self.widgetDescr['component'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':5, 'min':0, 'lockMin':1, 'lockMax':1, 'type':'int',
            'initialValue':0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'component'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='MSMSobject', name='MSMS')
        ip.append(datatype='TreeNodeSet',
                  balloon='Atoms for which triangles are extracted',
                  required=False, name='nodes')
        ip.append(datatype='int', required=False, name='selnum', defaultValue=3)
        ip.append(datatype='int', required=False, name='component', defaultValue=0)

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='vertices')
        op.append(datatype='faceIndices', name='indices')
        op.append(datatype='normals3D', name='vnormals')

        code = """def doit(self, msms, nodes, selnum, component):
    if self.inputPorts[3].widget:
        self.inputPorts[3].widget.configure(max=msms.rsr.nb-1, type='int')
    if nodes is None:
        if hasattr(msms, 'atoms'):
            atoms = msms.atoms
        else:
            atoms = None
            atomIndices = range(msms.nbat)
    else:
        atoms = nodes.findType(Atom)
    if selnum is None:
        selnum=3
    if component is None:
        component=0
    # get the index of every atom for which we want to display
    # the surface in the list of atoms used to compute the surface

    # special case all atoms
    atind=None
    if atoms and len(atoms)!=len(msms.atoms):
        inddict = {}
        i = 0
        for a in msms.atoms:
            inddict[a] = i
            i = i+1
        atind = map( lambda x, d=inddict: d[x], atoms)

    vf,vi,f = msms.getTriangles(atomIndices=atind, selnum=selnum,
                                component=component)
    normals = vf[:,3:6]
    faces = f[:,:3]
    if component > 0:
        normals = (-1.*Numeric.array( normals )).astype('f')
        faces = map( lambda x: (x[0], x[2], x[1]), f)

    self.outputData(vertices=vf[:,:3])
    self.outputData(indices=faces)
    self.outputData(vnormals=normals)\n"""

        self.setFunction(code)
        

    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)
        

class GetBuriedMSMS(NetworkNode):
    """Extract triangles from a given MSMS component that are burried by
a set of atoms from another molecule.

Input ports:
    MSMS: an MSMSobject instance as output by the MSMS calculation node
    nodes: a  TreeNodeSet describing a set of atoms. These atoms will be
            used to compute the part of the surface now rendered
            inaccessible to solvent.
    nodes2 :a TreeNodeSet describing parts of a molecule. This will be
           expanded to an AtomSet. Triangles for which at least 'selnum'
           vertices belong to atoms in this set will be extracted.
           If ommited, the whole surface is extracted.
    selnum: number of vertices per triangle that need to belong to atoms
            in nodes for a triangle to be output.
    component: MSMS component (i.e. a closed surface) to be extracted.
               Component 0 is always the external component.
    negate: When set to 1, the non-buried surface is extracted
    areaMode: calculation mode for surface area.
              Numeric: area computed using only flat triangles
              SemiAna: area computed using flat triangles for toroidal
                 patches and spherical triangles for contact and reentrant
                 patches.
OutputPort:
    vertices: array of X,Y,Z coordinates
    indices:  array of 0-based i,j,k indices describing triangular facets
    vnormals: array of Vx,Vy,Vz vertices normal vectors
    totalArea: surface area computed using only flat triangles. This port
               outputs a tuple of values for the solvent excluded and
               the solvent accessible surface areas (sesArea, sasAera).
    partialArea: surface area for the parts of the surface selected by nodes2
"""

    def __init__(self, name='get burried MSMS triangles', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['selnum'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':3, 'min':1, 'max':3, 'lockBMax': 1, 'lockMin': 1,
            'type':'int',
            'initialValue':3,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'nbVert'},
            }
        
        self.widgetDescr['component'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':5, 'min':0, 'lockMin':1, 'lockMax':1, 'type':'int',
            'initialValue':0,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'component'},
            }

        self.widgetDescr['negate'] = {
            'class':'NECheckButton',
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'negate'},
            }

        self.widgetDescr['areaMode'] = {
            'class':'NEComboBox',
            'choices':['Numeric', 'SemiAna'],
            'fixedChoices':True,
            'entryfield_entry_width':8,
            'labelGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'area mode'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='MSMSobject', name='MSMS')
        ip.append(datatype='TreeNodeSet',
                  balloon='nodes specifying atoms burying the surface',
                  name='nodes')
        ip.append(datatype='TreeNodeSet',
                  balloon='Atoms for which triangles are extracted',
                  required=False, name='nodes2')
        ip.append(datatype='int', required=False, name='selnum', defaultValue=3)
        ip.append(datatype='int', required=False, name='component', defaultValue=0)
        ip.append(datatype='int', required=False, name='negate', defaultValue=False)
        ip.append(datatype='string', required=False, name='areaMode', defaultValue='SemiAna')

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='vertices')
        op.append(datatype='faceIndices', name='indices')
        op.append(datatype='normals3D', name='vnormals')
        op.append(name='totalArea', datatype='tuple',
            balloon='outputs the total buried surface areas (sesArea,sasArea)')
        op.append(name='partialArea', datatype='tuple',
                  balloon='outputs the buried surface areas (sesArea,sasArea) for the atoms specfied in nodes2')

        code = """def doit(self, msms, nodes, nodes2, selnum, component, negate, areaMode):
    if not (msms or nodes):
        return
    self.inputPorts[3].widget.configure(max=msms.rsr.nb-1, type='int')
    if selnum is None:
        selnum=3
    if component is None:
        component=0
    # get the index of every atom for which we want to display
    # the surface in the list of atoms used to compute the surface
    atoms = nodes.findType(Atom)

    # identifyburied vertices
    msms.buriedVertices(atoms.coords, atoms.radius, component)

    # get the total buried surface area
    mode = _msms.MS_SEMI_ANALYTICAL
    if areaMode=='Numeric':
        mode = _msms.MS_NUMERICAL

    nareas = msms.buriedSurfaceArea(component, mode=mode)
    self.outputData(totalArea = (nareas['ses'][component],
                                 nareas['sas'][component]) )

    # get the triangles for the part of the buried surface corresponding
    # to the atomps specified by nodes2
    atind=None
    if nodes2 is not None:
        atoms2 = nodes2.findType(Atom)
        # special case all atoms
        if len(atoms2)!=len(msms.atoms):
            inddict = {}
            i = 0
            for a in msms.atoms:
                inddict[a] = i
                i = i+1
            atind = map( lambda x, d=inddict: d[x], atoms2)

    vf,vi,f = msms.getBuriedSurfaceTriangles(atomIndices=atind,
                  selnum=selnum, negate=negate, component=component)


    # vompute the partial surface area by summing up per vertex
    # buried surface areas

    # create a dict with buried vertex indices as keys
    usedvDict = {}
    for face in f:
        usedvDict[face[0]]=1
        usedvDict[face[1]]=1 
        usedvDict[face[2]]=1 

    # loop over all vertices and sum areas for the ones appearing in
    # the faces list
    psa1 = 0.0
    psa2 = 0.0
    for i in xrange(len(vf)):
         if usedvDict.has_key(i):
            psa1 = psa1 + vf[i][6]
            psa2 = psa2 + vf[i][7]

    self.outputData(partialArea = (psa1, psa2) )

    normals = vf[:,3:6]
    faces = f[:, :3]
    if component > 0:
        normals = (-1.*Numeric.array( normals )).astype('f')
        faces = map( lambda x: (x[0], x[2], x[1]), f)

    self.outputData(vertices=vf[:,:3])
    self.outputData(indices=faces)
    self.outputData(vnormals=normals)\n"""

        self.setFunction(code)
        

    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)


class AtomsProperty(NetworkNode):

    def __init__(self, name='Atoms Properties', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        #self.readOnly = 1

        self.widgetDescr['propertyName'] = {
            'class':'NEEntry', 'master':'node', 'width':13,
            'labelCfg':{'text':'prop. name:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atoms')
        ip.append(datatype='string', name='propertyName')

        op = self.outputPortsDescr
        op.append(datatype='list', name='propertyValues')

        code = """def doit(self, atoms, propertyName):
    import string
    if len(propertyName)==0: return
    values = []
    try:
        i = string.index(propertyName, '[')
        name = propertyName[:i]
        j = string.index(propertyName, ']')
        index = int(propertyName[i+1:j])
        values = eval('atoms.'+name)
        values = map( lambda x,i=index: x[i], values )
    except:
        try:
            values = eval('atoms.'+propertyName)
        except:
            import traceback
            traceback.print_stack()
            traceback.print_exc()
    if len(values):
        self.outputData(propertyValues=values)
"""

        self.setFunction(code)



class BondsGeometry(NetworkNode):

    def __init__(self, name='BondsGeometry', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atoms')

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='coords')
        op.append(datatype='faceIndices', name='indices')
        op.append(datatype='coordinates3D', name='nobonds')
        
        code = """def doit(self, atoms):
    coords = atoms.coords
    bnds, nobnds = atoms.bonds
    atoms._bndNb = range(len(atoms))
    bondIndices = map(lambda x: (x.atom1._bndNb,x.atom2._bndNb), bnds) 

    if coords:
        self.outputData(coords=coords)
    if bondIndices:
        self.outputData(indices=bondIndices)
    if nobnds:
        self.outputData(nobonds=nobnds)

    #clean up the atom set:
    for a in atoms:
        del a._bndNb
"""
            
        self.setFunction(code)


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)


class BondsByDist(NetworkNode):

    def __init__(self, name='buildBondsByDistance', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)

        self.widgetDescr['cut_off'] = {
            'class':'NEDial', 'master':'node', 'size':50,
            'oneTurn': 4.0,
            'initialValue':1.85,
            'labelCfg':{'text':'cutoff'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='MoleculeSet', name='molecules')
        ip.append(datatype='float', name='cut_off')

        op = self.outputPortsDescr
        op.append(datatype='MoleculeSet', name='molecules')

        code = """def doit(self, molecules, cut_off):
    if cut_off is None: cut_off=1.85

    for mol in molecules:
        mol.buildBondsByDistance(cut_off)
    self.outputData(molecules=molecules)\n"""

        if code: self.setFunction(code)

##################################################
# support ofr MSMS server node
#################################################

import types, numpy.oldnumeric as Numeric, re

class PublicHttpServer:

    def __init__(self, url):
        self.url = url
        self.arguments = {}


    def argCheck(self):
        return 1
    
    def run(self):
        import urllib
        if self.argCheck() == 0:
            return
        args = urllib.urlencode( self.arguments )
        f = urllib.urlopen(self.url, args)
        return f.read()


    def setArguments(self, **kw):
        from string import split
        for k,v in kw.items():
            if k not in self.arguments.keys():
                raise RuntimeError (k+' is not a valid argument')
            
            if type(v) in (types.FloatType, types.IntType,
                           types.LongType, types.StringType):
                self.arguments[k] = v
            else:
                pat = re.compile("[\[,\]\012]")
                arrayStr = re.sub(pat, '', str(Numeric.array(v).ravel()) )
                c = split(arrayStr)
                c = map( float, c )
                c = re.sub(pat, '', str(c) )
                self.arguments[k] = c
                
            
class MsmsServer(PublicHttpServer):

    def __init__(self, url):
        PublicHttpServer.__init__(self, url)
        self.arguments = {
            "Probe_Radius":1.5,
            "Density":1.0,
            "IsoRS":0,
            "Area":0,
            "All_Components":0,
            "R":0.0,
            "G":0.0,
            "B":0.0,
            "Coordinates": [],
            "Rendering":"Solid",
            "Lwidth":1.0,
            "Pwidth":1.0
            }

    def argCheck(self):
        return len(self.arguments['Coordinates'])


    def getTriangles(self, result):
        from string import split
        import re
        pat = re.compile(',')
        res = re.sub(pat, '', result)
        lines = split(res, '\n')
        coords = []
        triangles = []
        state = None
        for l in lines:
            if len(l)==0: continue
            if l[-1]=='}': state=None
            w = split(l[:-1])
            if len(w)==0: continue
            if w[0]=='point': state='coords'
            elif w[0]=='coordIndex': state='triangles'
            else:
                if state == 'coords':
                    coords.append( map(float, w) )
                elif state == 'triangles':
                    triangles.append( map(int, w[:3]) )
        return coords, triangles


class RemoteMSMS(NetworkNode):

    def __init__(self, name='Remote MSMS', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from DejaVu.IndexedPolygons import IndexedPolygons
        self.srf = IndexedPolygons('MSMS', inheritMaterial=0)

        self.widgetDescr['density'] = {
            'class':'NEDial', 'size':50, 'oneTurn':10, 'min':1.0,
            'increment':0.1, 'type':'float',
            'initialValue':1.0,
            'labelCfg':{'text':'density'},
            }
        
        self.widgetDescr['probe_radius'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':10, 'min':0.5, 'increment':0.1, 'type':'float',
            'initalValue':1.5,
            'labelCfg':{'text':'probe_radius'},
            }
        
        ip = self.inputPortsDescr
        ip.append(datatype='AtomSet', name='atoms')
        ip.append(datatype='vector', required=False, name='radii')
        ip.append(datatype='vector', required=False, name='colors')
        ip.append(datatype='float', required=False, name='density', defaultValue=1.)
        ip.append(datatype='float', required=False, name='probe_radius', defaultValue=1.5)

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='vertices')
        op.append(datatype='faceIndices', name='indices')
#        op.append(datatype='normals3D', name='vnormals')

        self.serverObj = MsmsServer(
            "http://www.scripps.edu/cgi-bin/sanner/msmsdemo.cgi")
        
        code = """def doit(self, atoms, radii, colors, density, probe_radius):
    if atoms is None: return        
    if radii is None:
        try: radii = atoms.radius
        except:
            print 'No Radii'
            return
    if colors is None:
        try: colors = atoms.colors
        except: colors = ( (0.7, 0.7, 0.7), )

    import numpy.oldnumeric as Numeric
    radii = Numeric.array(radii)
    radii.shape = (-1, 1)
    coords = Numeric.concatenate( (atoms.coords, radii), 1 )
    print coords.shape
    self.serverObj.setArguments( Coordinates=coords,
                                 Probe_Radius=probe_radius,
                                 Density=density )
    result = self.serverObj.run()
    v, f = self.serverObj.getTriangles(result)

    self.outputData(vertices=v)
    self.outputData(indices=f)\n"""

        self.setFunction(code)
        

class CrystalToCartesian(NetworkNode):
    """Converts crystall cell coordinates into cartesian coordinates.   """
    
    def __init__(self, name='Crystal To Cartesian', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        code = """def doit(self, coords, crystalCellCoords):
    if coords and crystalCellCoords:
        from mglutil.math.crystal import Crystal
        cc = crystalCellCoords
        self.crystal = Crystal(cc[0], cc[1])
        transformedCoords = self.crystal.toCartesian(coords)
        if transformedCoords is not None:
            self.outputData(transformedCoords=transformedCoords,
                            Crystal=self.crystal)
"""

        self.setFunction(code)

        ip = self.inputPortsDescr
        ip.append(datatype='coordinates3D', name='coords')
        ip.append(datatype='tuple', name='crystalCellCoords')

        op = self.outputPortsDescr
        op.append(datatype='coordinates3D', name='transformedCoords')
        op.append(datatype='CrystToReal', name='Crystal')


    def beforeAddingToNetwork(self, net):
        # import vizlib
        importVizLib(net)


from NetworkEditor.macros import MacroNode
class LinesMacro(MacroNode):
    def __init__(self, constrkw={}, name='lines macro', **kw):
        kw['name'] = name
        apply( MacroNode.__init__, (self,), kw)

    def beforeAddingToNetwork(self, net):
        MacroNode.beforeAddingToNetwork(self, net)
        # loading library molkitlib
        importMolKitLib(net)
        # loading library vizlib
        importVizLib(net)


    def afterAddingToNetwork(self):
        MacroNode.afterAddingToNetwork(self)
        node0 = self
        from DejaVu.VisionInterface.DejaVuNodes import vizlib
        node1 = node0.macroNetwork.ipNode
        node1.move(190, 10)
        node2 = node0.macroNetwork.opNode
        node2.move(210, 251)
        from MolKit.VisionInterface.MolKitNodes import BondsByDist
        node3 = BondsByDist(
            constrkw = {}, name='BondsByDist', library=molkitlib)
        node0.macroNetwork.addNode(node3,207,81)
        from MolKit.VisionInterface.MolKitNodes import NodeSelector
        node4 = NodeSelector(
            constrkw = {}, name='Select Nodes', library=molkitlib)
        node0.macroNetwork.addNode(node4,207,138)
        node4.inputPorts[2].widget.set("",0)
        from MolKit.VisionInterface.MolKitNodes import BondsGeometry
        node5 = BondsGeometry(
            constrkw = {}, name='BondsGeometry', library=molkitlib)
        node0.macroNetwork.addNode(node5,393,96)
        from DejaVu.VisionInterface.DejaVuNodes import IndexedPolylinesNE
        node6 = IndexedPolylinesNE(
            constrkw = {}, name='IndexedPolylines', library=vizlib)
        node0.macroNetwork.addNode(node6,393,186)

        ## saving connections for network Lines Macro ##
        node0.macroNetwork.connectNodes(
            node1, node3, "new", "molecules", blocking=True)
        node0.macroNetwork.connectNodes(
            node3, node4, "molecules", "nodes", blocking=True)
        node0.macroNetwork.connectNodes(
            node4, node5, "nodes", "atoms", blocking=True)
        node0.macroNetwork.connectNodes(
            node5, node6, "coords", "coords", blocking=True)
        node0.macroNetwork.connectNodes(
            node5, node6, "indices", "indices", blocking=True)
        node0.macroNetwork.connectNodes(
            node6, node2, "indexedPolylines", "new", blocking=True)

        ## modifying MacroOutputNode dynamic ports
        node2.inputPorts[1].configure(singleConnection=True)

        node0.shrink()
        ## reset modifications ##
        node0.resetTags()        
        node0.buildOriginalList()
        

from NetworkEditor.macros import MacroNode
class CPKMacro(MacroNode):
    def __init__(self, constrkw={}, name='CPK macro', **kw):
        kw['name'] = name
        apply( MacroNode.__init__, (self,), kw)

    def beforeAddingToNetwork(self, net):
        MacroNode.beforeAddingToNetwork(self, net)
        # loading library molkitlib
        importMolKitLib(net)
        # loading library vizlib
        importVizLib(net)
        # loading library stdlib
        from Vision.StandardNodes import stdlib
        net.editor.addLibraryInstance(stdlib, 'Vision.StandardNodes',
                                      'stdlib')

    def afterAddingToNetwork(self):
        MacroNode.afterAddingToNetwork(self)
        node0 = self
        from DejaVu.VisionInterface.DejaVuNodes import vizlib
        from MolKit.VisionInterface.MolKitNodes import AssignRadii
        node1 = node0.macroNetwork.ipNode
        node2 = node0.macroNetwork.opNode
        node3 = AssignRadii(
            constrkw = {}, name='Assign Radii', library=molkitlib)
        node0.macroNetwork.addNode(node3,76,83)
        node3.inputPorts[1].widget.set(1,0)
        from MolKit.VisionInterface.MolKitNodes import NodeSelector
        node4 = NodeSelector(
            constrkw = {}, name='Select Nodes', library=molkitlib)
        node0.macroNetwork.addNode(node4,78,213)
        node4.inputPorts[2].widget.set(".*",0)
        from MolKit.VisionInterface.MolKitNodes import AtomsProperty
        node5 = AtomsProperty(
            constrkw = {}, name='Extract Atom Property', library=molkitlib)
        node0.macroNetwork.addNode(node5,302,122)
        node5.inputPorts[1].widget.set("coords",0)
        from DejaVu.VisionInterface.DejaVuNodes import Spheres
        node6 = Spheres(constrkw = {}, name='spheres', library=vizlib)
        node0.macroNetwork.addNode(node6,211,225)
        node6.inputPortByName['radius'].unbindWidget()
        node6.inputPortByName['quality'].widget.set(10, run=False)

        ## saving connections for network CPK Macro ##
        node0.macroNetwork.connectNodes(
            node1, node3, "new", "molecules", blocking=True)
        node0.macroNetwork.connectNodes(
            node3, node4, "molecules", "nodes", blocking=True)
        node0.macroNetwork.connectNodes(
            node4, node5, "nodes", "atoms", blocking=True)
        node0.macroNetwork.connectNodes(
            node5, node6, "propertyValues", "coords", blocking=True)
        node0.macroNetwork.connectNodes(
            node6, node2, "spheres", "new", blocking=True)
        node0.macroNetwork.connectNodes(
            node3, node6, "radii", "radius", blocking=True)

        ## modifying MacroOutputNode dynamic ports
        node2.inputPorts[1].configure(singleConnection=True)

        node0.shrink()
        ## reset modifications ##
        node0.resetTags()
        node0.buildOriginalList()
        

from NetworkEditor.macros import MacroNode
class MSMSMacro(MacroNode):

    def __init__(self, constrkw={}, name='MSMS macro', **kw):
        kw['name'] = name
        apply( MacroNode.__init__, (self,), kw)

    def beforeAddingToNetwork(self, net):
        MacroNode.beforeAddingToNetwork(self, net)
        # loading library molkitlib
        importMolKitLib(net)
        # loading library vizlib
        importVizLib(net)

    def afterAddingToNetwork(self):
        MacroNode.afterAddingToNetwork(self)
        node0 = self
        from DejaVu.VisionInterface.DejaVuNodes import vizlib
        
        ## saving node MSMS Macro ##
        node1 = node0.macroNetwork.ipNode
        node2 = node0.macroNetwork.opNode
        from MolKit.VisionInterface.MolKitNodes import AssignRadii
        node3 = AssignRadii(
            constrkw = {}, name='Assign Radii', library=molkitlib)
        node0.macroNetwork.addNode(node3,91,84)
        from MolKit.VisionInterface.MolKitNodes import NodeSelector
        node4 = NodeSelector(
            constrkw = {}, name='Select Nodes', library=molkitlib)
        node0.macroNetwork.addNode(node4,94,149)
        node4.inputPorts[2].widget.set(".*",0)
        from MolKit.VisionInterface.MolKitNodes import AtomsAsMSMS
        node5 = AtomsAsMSMS(constrkw = {}, name='MSMS', library=molkitlib)
        node0.macroNetwork.addNode(node5,104,214)
        from MolKit.VisionInterface.MolKitNodes import GetMSMStriangles
        node6 = GetMSMStriangles(
            constrkw = {}, name='MSMS triang.', library=molkitlib)
        node0.macroNetwork.addNode(node6,373,102)
        from DejaVu.VisionInterface.DejaVuNodes import IndexedPolygonsNE
        node7 = IndexedPolygonsNE(
            constrkw = {}, name='indexedPolygons', library=vizlib)
        node0.macroNetwork.addNode(node7,383,171)
        ## saving connections for network MSMS Macro ##
        node0.macroNetwork.connectNodes(
            node3, node4, "molecules", "nodes", blocking=True)
        node0.macroNetwork.connectNodes(
                node4, node5, "nodes", "atoms", blocking=True)
        node0.macroNetwork.connectNodes(
                node5, node6, "MSMS", "MSMS", blocking=True)
        node0.macroNetwork.connectNodes(
                node6, node7, "vertices", "coords", blocking=True)
        node0.macroNetwork.connectNodes(
                node6, node7, "indices", "indices", blocking=True)
        node0.macroNetwork.connectNodes(
                node6, node7, "vnormals", "vnormals", blocking=True)
        node0.macroNetwork.connectNodes(
                node1, node3, "new", "molecules", blocking=True)
        node0.macroNetwork.connectNodes(
                node7, node2, "indexedPolygons", "new", blocking=True)

        ## modifying MacroOutputNode dynamic ports
        node2.inputPorts[1].configure(singleConnection=True)

        node0.shrink()
        ## reset modifications ##
        node0.resetTags()
        node0.buildOriginalList()
        

from Vision.VPE import NodeLibrary
molkitlib = NodeLibrary('MolKit', 'pink')

molkitlib.addNode(ReadMolecule, 'Read Molecule', 'Input')
molkitlib.addNode(AssignRadii, 'Assign Radii', 'Input')
molkitlib.addNode(AtomsProperty, 'Extract Atom Property', 'Input')
molkitlib.addNode(ReadMSMS, 'Read MSMS', 'Input')
molkitlib.addNode(PDB_SYMTRY, 'Get SYMTRY', 'Input')
molkitlib.addNode(PDB_MTRIX, 'Get MTRIX', 'Input')

molkitlib.addNode(NodeSelector, 'Select Nodes', 'Filter')
molkitlib.addNode(NodeSelector, 'Select MolFrag', 'Filter')

molkitlib.addNode(RMSDFromTwoAtomSet, 'RMSD', 'Mapper')
molkitlib.addNode(AtomsAsCPK, 'CPK', 'Mapper')
molkitlib.addNode(AtomsAsSticks, 'Sticks', 'Mapper')
molkitlib.addNode(BondsByDist, 'BondsByDist', 'Mapper')
molkitlib.addNode(BondsGeometry, 'BondsGeometry', 'Mapper')
molkitlib.addNode(CrystalToCartesian, 'Crystal To Cartesian', 'Mapper')

molkitlib.addNode(WritePDB, 'Write PDB', 'Output')
molkitlib.addNode(WriteMolecule, 'Write Molecule', 'Output')
molkitlib.addNode(writePDBFiles, 'Write PDB Files', 'Output')

molkitlib.addNode(LinesMacro, 'Lines Macro', 'Macro')
molkitlib.addNode(CPKMacro, 'CPK Macro', 'Macro')

molkitlib.addNode(CalculateGasteigerCharges, 'Calculate Gasteiger Charges', 'Mapper')
molkitlib.addNode(AddKollmanCharges, 'Add Kollman Charges', 'Mapper')
molkitlib.addNode(AddHydrogens, 'Add Hydrogens', 'Mapper')
molkitlib.addNode(AddPolarHydrogens, 'Add Polar Hydrogens', 'Mapper')
#    molkitlib.addNode(RemoteMSMS, 'remote MSMS', 'mapper')
#molkitlib.addNode(ParseCryst1, 'ParseCryst1', 'filter')
#molkitlib.addNode(CenterAndRadius, 'CenterAndRadius', 'filter')

if msmsFound==1:
    molkitlib.addNode(MSMSMacro, 'MSMS Macro', 'Macro')
    molkitlib.addNode(AtomsAsMSMS, 'MSMS', 'Mapper')
    molkitlib.addNode(SpheresAsMSMS, 'MSMS from spheres', 'Mapper')
    molkitlib.addNode(SaveMSMS, 'save MSMS', 'Output')
    molkitlib.addNode(GetMSMStriangles, 'MSMS triang.', 'Mapper')
    molkitlib.addNode(GetMSMSareas, 'MSMS areas', 'Mapper')
    molkitlib.addNode(GetBuriedMSMS, 'MSMS buried', 'Mapper')
    molkitlib.addNode(SESVertices, 'SESVertices', 'Filter')
    molkitlib.addNode(SurfaceAtoms, 'SurfaceAtoms', 'Filter')



UserLibBuild.addTypes(molkitlib, 'MolKit.VisionInterface.MolKitTypes')

try:
    UserLibBuild.addTypes(molkitlib, 'DejaVu.VisionInterface.DejaVuTypes')
except:
    pass
