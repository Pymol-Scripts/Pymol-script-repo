## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#!/usr/bin/env python2.3
#
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/MoleculePreparation.py,v 1.61.2.5 2009/09/25 15:45:27 rhuey Exp $ 
#
#
#
import os, string
import numpy.oldnumeric as Numeric, math 

from MolKit import Read
from MolKit.torTree import TorTree
from MolKit.pdbWriter import PdbqWriter, PdbqsWriter, PdbqtWriter
import MolKit.molecule
import MolKit.protein
from MolKit.protein import ResidueSet
from MolKit.mol2Parser import Mol2Parser

from MolKit.molecule import AtomSet, Bond, BondSet
from MolKit.chargeCalculator import GasteigerChargeCalculator 
from MolKit.chargeCalculator import KollmanChargeCalculator
from MolKit.hydrogenBuilder import HydrogenBuilder
#import the Kollman charges dictionary to get keys it recognizes...
from MolKit.data.qkollua import q
from MolKit.bondSelector import RotatableBondSelector, AmideBondSelector, LeafBondSelector
from MolKit.bondSelector import GuanidiniumBondSelector



from AutoDockTools.atomTypeTools import NonpolarHydrogenMerger
from AutoDockTools.atomTypeTools import LonepairMerger
from AutoDockTools.atomTypeTools import SolvationParameterizer
from AutoDockTools.atomTypeTools import AromaticCarbonManager
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper

from AutoDockTools.AutoDockBondClassifier import AutoDockBondClassifier



class AutoDockMoleculePreparation:
    """
    Base class for facades for preparing MolKit molecules to be used in an
    AutoDock experiment.

    This preparation involves:
        1. optional cleanup
        2. making atoms conform to AD atom format:
            adding charges
            merging nonpolarHydrogens (can be overridden)
            merging lonepairsHydrogens (can be overridden)
            adding solvation parameters
          [optional:
            adding Hydrogens]
        3. writing an pdbq[s,t] outputfile 


    Collaborators:
        in MolKit:
            for adding charges:
                KollmanChargeCalculator
                [optional: GasteigerChargeCalculator]
            for adding hydrogens:
                HydrogenBuilder
        in AutoDockTools:
            atomTypeManager classes for conforming to AD atomtypes:
                LonepairMerger
                NonpolarHydrogenMerger
        

    Preparation of a molecule 
        buildBonds if necessary
        adds 'autodock_element field' to all Atoms
        detects whether molecule isPeptide or not (for charge type to use)


    Methods to be shared by derived classes
        calcChargeError: calculation of totalCharge error
        detectIsPeptide: returns boolean isPeptide
        findNearest: calculation of atom nearest another (non-bonded) atom
"""

    
    std_types = ['CYS', 'ILE', 'SER', 'VAL', 'GLN', 'LYS', 'ASN', 'PRO', 'THR', 'PHE', 'ALA', 'HIS', 'GLY', 'ASP', 'LEU', 'ARG', 'TRP', 'GLU', 'TYR', 'MET', 'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']
    #what about nucleic acids???

    def __init__(self, molecule, mode='', repairs='',
                    charges_to_add=None, cleanup='',
                    outputfilename=None, debug=False,
                    version=3, delete_single_nonstd_residues=False):
        #debug = True
        self.debug = debug
        self.molecule = mol = molecule
        self.mode = mode
        self.lenNPHS = 0
        #fix atom elements if necessary:
        self.autodock_element_dict = {'Cl':'c', 'Fe': 'f', 'Br':'b'}
        self.reverse_autodock_elements = {'c':'Cl', 'f':'Fe','b':'Br'}
        self.version = version
        self.delete_single_nonstd_residues = delete_single_nonstd_residues
        if version==3:
            revlist = ['c','f','b']
            reversedAts = AtomSet(filter(lambda x: x.element in revlist, mol.allAtoms))
            for at in reversedAts:
                if debug: print 'restoring ', at.element, ' to ',
                at.autodock_element = at.element
                at.element = self.reverse_autodock_elements[at.element]
                if debug: print  at.element, ' for processing...'
        # build bonds if necessary      #ALWAYS DO THIS!!!
        if not len(mol.allAtoms.bonds[0]): #what if not enough bonds
            mol.buildBondsByDistance()

        #PROCESS molecule
        self.cleanup_type_list = string.split(cleanup, '_')
        #if mode=='automatic':  #DO this ANYWAY???
        #remove waters + chains composed of nonstd residues
        self.cleanUpResidues(mol, self.cleanup_type_list)

        # REPAIR: possibly add bonds to non-bonded atoms + add hydrogens
        self.repair_type_list = string.split(repairs, '_')
        #if mode=='automatic':  #SHOULD THIS BE DONE ANYWAY???
        self.repairMol(mol, self.repair_type_list)
        ##3/25/2008: have to remove lps first because there is no atomic
        #number for Lp or LP in PyBabel/babelElements:
        merge_lonepairs = 'lps' in self.cleanup_type_list
        if merge_lonepairs:
            LPM = self.LPM = LonepairMerger()
            #lps = filter(lambda x: x.element=='Xx' and \
                  #(x.name[0]=='L' or x.name[1]=='L'), mol.allAtoms)
            lenLPS = self.lenLPS = LPM.mergeLPS(mol.allAtoms)
            if debug: print 'merged ', lenLPS, ' lonepairs'
            self.cleanup_type_list.remove('lps')
        # CHARGES: ALWAYS add charges last after atoms are set but
        # before merging nphs and lps
        self.chargeType = charges_to_add
        self.chargeError = 'ERROR'
        if debug: print "charges_to_add=", charges_to_add
        #if mode=='automatic':  #SHOULD THIS BE DONE ANYWAY???
        self.addCharges(mol, charges_to_add)
        if self.chargeType!='ERROR':
            self.chargeError = self.calcChargeError()
        #need a two step cleanup to do merges after adding bonds/hydrogens
        #and AFTER adding charges
        # CLEANUP: possibly merge nonpolarhydrogens and/or lonepairs
        self.cleanUpAtomTypes(mol, self.cleanup_type_list)
        #set autodock_element for all atoms
        mol.allAtoms.autodock_element = mol.allAtoms.element
        self.outputfilename = outputfilename
        if debug: print "end of base class init"


    def repairMol(self, mol, repair_list):
        #@@ what about metals?
        # possibly check for non-bonded atoms and add bonds to them
        # possibly add hydrogens 
        if 'bonds' in repair_list:
            non_bonded_atoms = mol.allAtoms.get(lambda x: len(x.bonds)==0)
            if non_bonded_atoms:
                if self.debug: print "!!adding bonds to ", non_bonded_atoms.name, "!!"
                for a in non_bonded_atoms:
                    nearestAt = self.findNearest(a, mol.allAtoms)
                    Bond(a, nearestAt, bondOrder=1, origin='AddedBond', check=0,
                                addBond=1)
                    if self.debug: print 'added bond:', a.bonds[0]
            repair_list.remove('bonds')
        # possibly add hydrogens 
        #MODE SWITCH 1: adding hydrogens  ||IN USE||
        self.newHs = 0
        hs = mol.allAtoms.get(lambda x: x.element=='H')
        if 'hydrogens' in repair_list:
            if self.debug: print "addinghydrogens!"
            self.newHs = self.addHydrogens(mol)
            if self.debug: print "added ", self.newHs, " hydrogens"
            repair_list.remove('hydrogens')
        elif not hs and 'checkhydrogens' in self.repair_type_list:
            if self.debug: print "addinghydrogens from check hydrogens"
            self.newHs = self.addHydrogens(mol)
            if self.debug: print "added ", self.newHs, " hydrogens"
            repair_list.remove('checkhydrogens')

    def addHydrogens(self, mol):
        beforeLen = len(mol.allAtoms)
        #NB: this adds all hydrogens (?!?)
        HB = HydrogenBuilder()
        HB.addHydrogens(mol)
        afterLen = len(mol.allAtoms)
        newHs = afterLen - beforeLen
        #NB this adds the HYDROGENS as "ATOM" not HETATM
        return newHs


    def addCharges(self, mol, charges_to_add):
        #detect whether isPeptide first
        debug = self.debug
        if debug: print "in addCharges:mol.name=", mol.name, " and charges_to_add=", charges_to_add
        isPeptide = mol.isPeptide = self.detectIsPeptide()
        if charges_to_add=='gasteiger':
            chargeCalculator = self.chargeCalculator = GasteigerChargeCalculator()
            if isPeptide:
                print "adding gasteiger charges to peptide"
                chargeCalculator = self.chargeCalculator = GasteigerChargeCalculator()
        elif charges_to_add=='Kollman':
            if not mol.isPeptide:
                print "adding Kollman charges to non-peptide"
            chargeCalculator = self.chargeCalculator = KollmanChargeCalculator()
        #KEEP charges from file, if there are any
        #MODE SWITCH 2: adding charges????  ||NOT IN USE||
        if debug: print "self.chargeType=", self.chargeType
        if self.chargeType is None:
            len_zero_charges = len(filter(lambda x: x.charge==0, 
                                mol.chains.residues.atoms))
            if len_zero_charges==len(mol.allAtoms):
                print "WARNING: all atoms in '%s' had zero charges! Adding gasteiger charges..."%mol.name
                chargeCalculator = self.chargeCalculator = GasteigerChargeCalculator()
                chargeCalculator.addCharges(mol.allAtoms)
            else:
                len_charged = len(filter(lambda x:x.chargeSet!=None, 
                                mol.chains.residues.atoms))
                if len_charged!=len(mol.chains.residues.atoms):
                    #charges could come from file or ?
                    print "WARNING: some atoms in '%s' had no charge! Adding gasteiger charges to all..."%mol.name
                if debug: print "no change to charges!"
            self.chargeType = mol.allAtoms[0].chargeSet
        else: 
            chargeCalculator.addCharges(mol.allAtoms)
            if debug: print 'added ' + charges_to_add + 'charges to ', mol.name


    def cleanUpAtomTypes(self, mol, cleanup_list):
        debug = self.debug
        merge_nonpolar_hydrogens = 'nphs' in cleanup_list
        if merge_nonpolar_hydrogens:
            NPHM = self.NPHM = NonpolarHydrogenMerger()
            lenNPHS = self.lenNPHS = NPHM.mergeNPHS(mol.allAtoms)
            if debug: print 'merged ', lenNPHS, ' nonpolar hydrogens'
            cleanup_list.remove('nphs')
        merge_lonepairs = 'lps' in cleanup_list
        if merge_lonepairs:
            LPM = self.LPM = LonepairMerger()
            #lps = filter(lambda x: x.element=='Xx' and \
                  #(x.name[0]=='L' or x.name[1]=='L'), mol.allAtoms)
            lenLPS = self.lenLPS = LPM.mergeLPS(mol.allAtoms)
            if debug: print 'merged ', lenLPS, ' lonepairs'
            cleanup_list.remove('lps')
        if self.delete_single_nonstd_residues:
            #delete non-standard residues from within a single chain
            non_std = []
            for res in mol.chains.residues:
                if res.type not in self.std_types:
                    non_std.append(res)
            #delete non_std
            names = ""
            #remove all bonds to other atoms:
            for r in non_std:
                for a in r.atoms:
                    for b in a.bonds:
                        at2 = b.atom1
                        if at2==a: at2 = b.atom2
                        at2.bonds.remove(b)
                        a.bonds.remove(b)
                names = names + r.parent.id +r.name + '_'
            print "\'Deleting non-standard residues:" + names + " from " + mol.name
            for res in non_std:
                res.parent.remove(res)
                del res
            #fix allAtoms short cut
            mol.allAtoms = mol.chains.residues.atoms
        delete_alternateB = 'deleteAltB' in cleanup_list        
        if delete_alternateB:
            if self.debug: print "deleting ", len(mol.allAtoms), " alternateB positions"
            altB_ats = mol.allAtoms.get('*@B')
            if not len(altB_ats):
                return
            for a in altB_ats:
                res = a.parent
                for b in a.bonds:
                    at2 = b.atom1
                    if at2==a: at2 = b.atom2
                    at2.bonds.remove(b)
                    a.bonds.remove(b)
                res.remove(a)
                del a
            #remove any remaining @A
            altA_ats = mol.allAtoms.get('*@A')
            for a in altA_ats:
                a.name = a.name[:-2]
            mol.allAtoms = mol.chains.residues.atoms
            if self.debug: print "resetting numbers to 1-", len(mol.allAtoms)+1
            mol.allAtoms.number = range(1, len(mol.allAtoms)+1)



    def cleanUpResidues(self, mol, cleanup_list):
        self.lenHOHS = 0
        remove_waters = 'waters' in cleanup_list
        if remove_waters:
            hohs = mol.allAtoms.get(lambda x: x.parent.type=='HOH')
            #this line added at Alex Perryman's suggestion:
            hohs2 = mol.allAtoms.get(lambda x: x.parent.type=='WAT')
            hohs = hohs + hohs2 
            if hohs:
                #remove(hohs)
                lenHOHS = self.lenHOHS = len(hohs)
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
                        if self.debug: print "removing residue: ", res.name
                        chain.remove(res)
                        if len(chain.children)==0:
                            mol = chain.parent
                            if self.debug: print "removing chain", chain.id
                            mol.remove(chain)
                            del chain
                        del res
                    del h 
                #fix allAtoms short cut
                mol.allAtoms = mol.chains.residues.atoms    
            cleanup_list.remove('waters')

        remove_nonstdres = 'nonstdres' in cleanup_list
        self.lenChainsRemoved = 0
        if remove_nonstdres: 
            if len(mol.chains)>1:
                chains_to_delete = []
                for c in mol.chains:
                    num_res = len(c.residues)
                    non_std = []
                    for res in c.residues:
                        if res.type not in self.std_types:
                            non_std.append(res)
                    if len(c.residues)==len(non_std):
                        print "\'"+ c.name, "\' apparently composed of not std residues. Deleting "
                        chains_to_delete.append(c)
                if len(chains_to_delete):
                    self.lenChainsRemoved = len(chains_to_delete)
                    for c in chains_to_delete:
                        #remove(c)     #NB this would remove HOH and ligand from 1HSG
                        c.parent.remove(c)
                        del c
                    #fix allAtoms short cut
                    mol.allAtoms = mol.chains.residues.atoms    
            #if self.delete_single_nonstd_residues:
            #    #delete non-standard residues from within a single chain
            #    non_std = []
            #    for res in mol.chains.residues:
            #        if res.type not in self.std_types:
            #            non_std.append(res)
            #    #delete non_std
            #    names = ""
            #    #remove all bonds to other atoms:
            #    for r in non_std:
            #        for a in r.atoms:
            #            for b in a.bonds:
            #                at2 = b.atom1
            #                if at2==a: at2 = b.atom2
            #                at2.bonds.remove(b)
            #                a.bonds.remove(b)
            #        names = names + r.parent.id +r.name + '_'
            #    print "\'Deleting non-standard residues:" + names + " from " + mol.name
            #    for res in non_std:
            #        res.parent.remove(res)
            #        del res
            #    #fix allAtoms short cut
            #    mol.allAtoms = mol.chains.residues.atoms
            cleanup_list.remove('nonstdres')


    def calcChargeError(self):
        s = Numeric.add.reduce(self.molecule.allAtoms.charge)
        return min(math.ceil(s)-s, s-math.floor(s))


    def detectIsPeptide(self):
        isPeptide = True
        d = {}
        for r in self.molecule.chains.residues:
            d[r.type] = 0
        for t in d.keys():
            if t not in q.keys():
                isPeptide = False
                break
        return isPeptide


    def findNearest(self, atom, bonded_atoms):
        lenK = len(bonded_atoms)
        lenC = 1
        nonbonded_atoms = AtomSet([atom])
        c = Numeric.array(nonbonded_atoms.coords, 'f')
        k = Numeric.array(bonded_atoms.coords, 'f')
        bigK = Numeric.resize(k, (lenK, lenK, 3))
        c.shape = (1, 1, 3)
        bigM = bigK[:1]
        d = bigM - c
        dSQ = d*d
        dSQMAT = Numeric.sum(dSQ, 2)
        mm = dSQMAT[0][0]
        mindex = 0
        for i in range(lenK):
            if dSQMAT[0][i]<mm:
                mm = dSQMAT[0][i]
                mindex = i
        #print "found closest atom %d at dist %f" %( mindex, mm)
        return bonded_atoms[mindex]



class ReceptorPreparation(AutoDockMoleculePreparation):
    """
    Facade for preparing a MolKit molecule to be used as the receptor in an
    AutoDock experiment derived from AutoDockMoleculePreparation class.


    Receptor preparation involves:
        1. adding solvation parameters
        2. writing a pdbqs outputfile 

    Receptor preparation Collaborators extend base class collaborators:
        in AutoDockTools:
            in atomTypeManager classes for conforming to AD atomtypes:
                SolvationParameterizer
       in this file:
            ReceptorWriter
        
"""

    def __init__(self, molecule, mode='automatic', repairs='checkhydrogens',
                    charges_to_add='Kollman',
                    cleanup='nphs_lps_waters_nonstdres',
                    outputfilename=None, debug=False,
                    preserved={},
                    delete_single_nonstd_residues=False):

        AutoDockMoleculePreparation.__init__(self, molecule, mode, repairs,
                                            charges_to_add, cleanup,
                                            debug=debug,
                                            delete_single_nonstd_residues=False)

        if len(preserved):
            for atom, chargeList in preserved.items():
                atom._charges[chargeList[0]] = chargeList[1]
                atom.chargeSet = chargeList[0]
                if debug: print 'preserved charge on ', atom.name, ' charge=', atom.charge

        # add solvation parameters here:
        SP = self.SP = SolvationParameterizer()
        unknownatoms = SP.addParameters(molecule.chains.residues.atoms)
        if len(unknownatoms) and debug:
            print len(unknownatoms), " unknown atoms: set solvation parameters to 0. 0."

        self.writer = ReceptorWriter()
        #MODE SWITCH 5: write outputfile     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: what records should be written?
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def summarize(self):
        mol = self.molecule
        msg = "setup "+ mol.name+ ":\n"
        if self.newHs!=0:   
            msg = msg + " added %d new hydrogen(s)\n" %self.newHs
        if self.chargeType!=None:   
            if self.chargeType in ['pdbq', 'mol2', 'pdbqt', 'pdbqs']:
                msg = msg + " kept charges from " + self.chargeType + " file\n"
            else:
                msg = msg + " added %s charges\n" %self.chargeType
        if self.lenNPHS:
            msg = msg + " merged %d non-polar hydrogens\n" %self.lenNPHS
        if self.lenLPS:
            msg = msg + " merged %d lone pairs\n" %self.lenLPS
        totalCh = Numeric.add.reduce(mol.allAtoms.charge)
        if hasattr(self, 'chargeError') and abs(self.chargeError)>0.000005:
            msg = msg + " total charge error = %6.4f\n"%self.chargeError
        return msg


    def write(self, outputfilename=None):
        #write to outputfile
        mol = self.molecule
        if not outputfilename and not self.outputfilename:
            outputfilename = self.outputfilename = mol.name + ".pdbqs"
            if self.debug: print 'using std outputfilename ', outputfilename
        #should this be done at Molecule level or at Atom level?
        #Atom level would allow multiple molecules in Receptor, at no price
        self.writer.write(mol, outputfilename)
        if self.debug: print "wrote ", mol.name, " to ", outputfilename



class AD4ReceptorPreparation(AutoDockMoleculePreparation):
    """
    Facade for preparing a MolKit molecule to be used as the receptor in an
    AutoDock4 experiment derived from AutoDockMoleculePreparation class.

    Receptor preparation involves:
        1. adding autodock_element types
        2. writing a pdbqt outputfile 

    Receptor preparation Collaborators extend base class collaborators:
        in AutoDockTools:
            in atomTypeManager classes for conforming to AD atomtypes:
                AutoDock4_AtomTyper
       in this file:
            AD4ReceptorWriter
        
"""

    def __init__(self, molecule, mode='automatic', repairs='checkhydrogens',
                    charges_to_add='gasteiger',
                    cleanup='nphs_lps_waters_nonstdres',
                    outputfilename=None, debug=False,
                    version=4, preserved={},
                    delete_single_nonstd_residues=False):


        AutoDockMoleculePreparation.__init__(self, molecule, 
                        mode=mode, repairs=repairs, 
                        charges_to_add=charges_to_add, cleanup=cleanup,
                        outputfilename=outputfilename, debug=debug,
                        version=version, delete_single_nonstd_residues=delete_single_nonstd_residues)

        #NEED TO TYPE ATOMS!!
        try:
            delattr(mol.allAtoms, 'autodock_element')
        except:
            pass
        atomTyper = AutoDock4_AtomTyper()
        atomTyper.setAutoDockElements(molecule, reassign=True) #4/15/05 catch here?
        if len(preserved):
            for atom, chargeList in preserved.items():
                atom._charges[chargeList[0]] = chargeList[1]
                atom.chargeSet = chargeList[0]
                if debug: print 'preserved charge on ', atom.name, ' charge=', atom.charge

        self.writer = AD4ReceptorWriter()
        #MODE SWITCH 5: write outputfile     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: what records should be written?
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def write(self, outputfilename=None):
        #write to outputfile
        mol = self.molecule
        if not outputfilename and not self.outputfilename:
            outputfilename = self.outputfilename = mol.name + ".pdbqt"
            if self.debug: print 'using std outputfilename ', outputfilename
        #should this be done at Molecule level or at Atom level?
        #Atom level would allow multiple molecules in Receptor, at no price
        self.writer.write(mol, outputfilename)
        if self.debug: print "wrote ", mol.name, " to ", outputfilename

    


class ReceptorWriter:
    """
    ReceptorWriter class writes a receptor molecule which has a ReceptorPreparationObject to a file
"""

    def __init__(self, write_CONECT=False):
        self.writer = PdbqsWriter()
        #each atom must have charge and AtVol and AtSolPar
        self.conditions = [lambda x: hasattr(x, 'AtVol'), \
                           lambda x: hasattr(x, 'AtSolPar'), \
                           lambda x: x.chargeSet is not None]
        self.write_CONECT = write_CONECT


    def write(self, receptor, outputfile):
        #should this be done at Molecule level or at Atom level?
        receptor_atoms = receptor.chains.residues.atoms
        len_receptor_atoms = len(receptor_atoms)

        ##check each condition for receptor_atoms
        for cond in self.conditions:
            assert len(filter(cond, receptor_atoms))==len_receptor_atoms

        outptr = open(outputfile, 'w')
        records = ['ATOM']
        #if self.write_CONECT:
        #    records.append('CONECT')
        #for at in receptor_atoms:
        #    self.writer.write_atom(outptr, at)
        for ch in receptor.chains:
            for at in ch.residues.atoms:
                self.writer.write_atom(outptr, at)
            rec = self.writer.defineTERRecord(at)
            outptr.write(rec)
        #put in the TER and END lines???
        #optional CONECT records???
        if self.write_CONECT:
            self.writeCONECTRecords(receptor, outptr)
        outptr.close()
        

    def writeCONECTRecords(self, fptr, mol):
        atms = mol.allAtoms
        for atm in atms:
            rec = 'CONECT%5i'%atm.number
            for b in atm.bonds:
                a2 = b.atom1
                if a2==atm: a2 = b.atom2
                if not a2 in atms: continue #don't write intermolecular bonds
                rec = rec + '%5i'%a2.number 
            rec = rec + '\n'
            fptr.write(rec)




class AD4ReceptorWriter(ReceptorWriter):
    """
    AD4ReceptorWriter class writes a receptor molecule which has a ReceptorPreparationObject to a pdbqt file
"""

    def __init__(self):
        ReceptorWriter.__init__(self)
        self.writer = PdbqtWriter()
        self.conditions = [lambda x: hasattr(x, 'autodock_element'), \
                           lambda x: x.chargeSet is not None]



class LigandPreparation(AutoDockMoleculePreparation):

    """

    Facade for preparing a MolKit molecule to be used as the ligand in an
    AutoDock experiment.

    LigandPreparation involves:
        1. managing definition of flexibility pattern in ligand
           using metaphor of a Tree
            +setting ROOT
            +setting pattern of rotatable bonds
            +setting TORSDOF, number of possible rotatable bonds-
                          number of bonds which only rotate hydrogens
           optional:
            -disallowing amide torsions
            -disallowing peptidebackbone torsions
            -disallowing guanidinium torsions
            +/- toggling activity of any rotatable bond
        2. writing an outputfile with keywords recognized by AutoDock
        3. writing types list and number of active torsions to a dictionary file

    LigandPreparation Collaborators:
        in AutoDockTools:
            for rotability pattern management:
                RotatableBondManager
                [which uses an AutoDockTools/AutoDockBondClassifier
                for detecting possible flexibility patterns in bonds] 
            for managing cyclic aromatic carbons:
                AromaticCarbonManager
                (including changing cutoff from 7.5)
        
"""

    def __init__(self, molecule, mode='automatic', repairs='hydrogens_bonds',
                    charges_to_add='gasteiger', cleanup='nphs_lps',
                    allowed_bonds='backbone',  #backbone on, amide + guanidinium off
                    root='auto', outputfilename=None, dict=None, debug=False,
                    check_for_fragments=False, bonds_to_inactivate=[],
                    inactivate_all_torsions=False,
                    version=3, limit_torsions=False,
                    delete_single_nonstd_residues=False,
                    detect_bonds_between_cycles=False):
        # why aren't repairs, cleanup and allowed_bonds lists??
        # to run tests: use allowed_bonds = 'guanidinium'
        if debug: print "LPO: allowed_bonds=", allowed_bonds
        self.detect_bonds_between_cycles = detect_bonds_between_cycles


        AutoDockMoleculePreparation.__init__(self, molecule, mode, repairs,
                                            charges_to_add, cleanup,
                                            debug=debug, version=version,
                                            delete_single_nonstd_residues=False)

        self.version = version
        self.charge_error_tolerance = 0.000005
        #FOR FLEXIBILITY MODEL: setup RotatableBondManager
        #process bonds
        #MODE SWITCH 4: disallow allowed_bonds ||NOT IN USE||
        self.RBM = RotatableBondManager(molecule,
                            allowed_bonds, root, debug=False,
                            check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate,
                            detectAll=self.detect_bonds_between_cycles)
        if inactivate_all_torsions is True:
            #in this case, make all active torsions inactive
            self.RBM.set_all_torsions(False)
        elif limit_torsions is not False:
            self.RBM.limit_torsions(limit_torsions)

        #detect aromatic cycle carbons and rename them for AutoDock3 
        rename = self.version==3
        ACM = self.ACM = AromaticCarbonManager(rename=rename)
        self.aromCs = self.ACM.setAromaticCarbons(molecule) 

        #optional output summary filename to write types and number of torsions
        self.dict = dict

        self.outputfilename = outputfilename
        #if mode is 'automatic': write outputfile now 
        #without waiting to set torsions etc interactively
        #MODE SWITCH 5: write outputfile     ||IN USE||
        self.writer = LigandWriter()
        if mode=='automatic':
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def summarize(self):
        mol = self.molecule
        msg = "setup "+ mol.name+ ":\n"
        if self.newHs!=0:   
            msg = msg + " added %d new hydrogen(s)\n" %self.newHs
        if self.chargeType!=None:   
            if self.chargeType in ['pdbq', 'mol2']:
                msg = msg + " kept charges from " + self.chargeType + " file\n"
            else:
                msg = msg + " added %s charges\n" %self.chargeType
        if self.lenNPHS:
            msg = msg + " merged %d non-polar hydrogens\n" %self.lenNPHS
        if self.lenLPS:
            msg = msg + " merged %d lone pairs\n" %self.lenLPS
        if len(self.aromCs):
            msg = msg + " found %d aromatic carbons\n" %len(self.aromCs)
        totalCh = Numeric.add.reduce(mol.allAtoms.charge)
        msg = msg + " detected %d rotatable bonds\n" %len(mol.possible_tors_bnds)
        msg = msg + " set TORSDOF to %d\n"%mol.TORSDOF
        #if hasattr(self, 'chargeError') and abs(self.chargeError)>0.000005:
        if hasattr(self, 'chargeError') and abs(self.chargeError)>self.charge_error_tolerance:
            msg = msg + " total charge error = %6.4f\n"%self.chargeError
        return msg

    def set_autodock_element_types_for_writing(self, allAtoms):
        sublist = self.autodock_element_dict.keys()
        problem_ats = AtomSet(filter(lambda x: x.element in sublist, allAtoms))
        for at in problem_ats:
            #there is a substitute 'name' for this atom
            if self.debug: print "changed ", at.name,
            newval = self.autodock_element_dict[at.element]
            at.name = newval
            #have to change element to get correct written output
            #because the stupid pdbWriter writes the element not the name
            at.element = newval
            at.autodock_element = at.element
            if self.debug: print "to ", at.name
        return problem_ats


    def restore_autodock_element_types_after_writing(self, problem_ats):
        for at in problem_ats:
            #there is a substitute 'name' for this atom
            if self.debug: print "changed back", at.name,"'s element "
            newval = self.reverse_autodock_elements[at.element]
            #have to change element to get correct written output
            #because the stupid pdbWriter writes the element not the name
            at.element = newval
            if self.debug: print "to ", at.element


    def write(self, outputfilename=None):
        #write ligand to outputfile
        mol = self.molecule
        if not hasattr(mol, 'ROOT'):
            print 'must specify ROOT before writing'
            return 'ERROR' 
        if hasattr(mol, 'torTree') and not hasattr(mol, 'torscount'):
            mol.torscount = len(mol.torTree.torsionMap)
        if not outputfilename:
            stem = os.path.splitext(os.path.basename(mol.parser.filename))[0]
            if self.version>=4:
                outputfilename = stem + ".pdbqt"
            else:
                outputfilename = stem + ".pdbq"
            if self.debug: print 'using std outputfilename ', outputfilename
            #if self.debug: print 'using std outputfilename ', outputfilename
        #writer = LigandWriter()
        #print "calling writer.write with ", mol.name, ' ',  outputfilename
        #change the names of Cl, Fe and Br atoms
        problem_ats = self.set_autodock_element_types_for_writing(mol.allAtoms)
        #if self.version==3:
            #sublist = self.autodock_element_dict.keys()
            #problem_ats = AtomSet(filter(lambda x: x.element in sublist, mol.allAtoms))
            #for at in problem_ats:
                ##there is a substitute 'name' for this atom
                #if self.debug: print "changed ", at.name,
                #newval = self.autodock_element_dict[at.element]
                #at.name = newval
                ##have to change element to get correct written output
                ##because the stupid pdbWriter writes the element not the name
                #at.element = newval
                #at.autodock_element = at.element
                #if self.debug: print "to ", at.name
        self.writer.write(mol, outputfilename)
        self.outputfilename = outputfilename
        #print 'back from writing'
        if self.dict is not None:
            if not os.path.exists(self.dict):
                fptr = open(self.dict, 'a')
                ostr = "summary = d = {}\n"
                fptr.write(ostr)
            else:
                fptr = open(self.dict, 'a')
            type_dict = {}
            for a in self.molecule.allAtoms:
                type_dict[a.autodock_element] = 1
            atom_types = type_dict.keys()
            atom_types.sort()
            ostr = "d['" + self.molecule.name +"'] = {"
            ostr = ostr + "'atom_types': [" 
            for t in atom_types[:-1]:
                ostr = ostr + "'%s', "%t
            ostr = ostr + "'%s' "%atom_types[-1]
            ostr = ostr + "],\n\t\t\t'rbonds':" + str(self.molecule.torscount)
            ostr = ostr + ",\n\t\t\t'zero_charge' : ["
            zc_ats = self.molecule.allAtoms.get(lambda x: x.charge==0)
            for at in zc_ats:
                ostr = ostr + "'%s', "%at.name
            ostr = ostr + "],\n\t\t\t}\n"
            fptr.write(ostr)
            fptr.close()
        self.restore_autodock_element_types_after_writing(problem_ats)
        #if self.version==3:
            #for at in problem_ats:
                ##there is a substitute 'name' for this atom
                #if self.debug: print "changed back", at.name,"'s element "
                #newval = self.reverse_autodock_elements[at.element]
                ##have to change element to get correct written output
                ##because the stupid pdbWriter writes the element not the name
                #at.element = newval
                #if self.debug: print "to ", at.element


    def autoroot(self):
        self.RBM.autoroot()


    def setroot(self, index):
        self.RBM.setroot(index)


    def limit_torsions(self, val, type):
        self.RBM.limit_torsions(val, type)


    def set_all_torsions(self, flag):
        mol = self.molecule
        if self.debug: print "LPO:set_all_torsions with flag=", flag
        if self.debug: print "activeBonds=", len(filter(lambda x: x.activeTors, mol.allAtoms.bonds[0]))
        self.RBM.set_all_torsions(flag)
        if self.debug: print "2: activeBonds=", len(filter(lambda x: x.activeTors, mol.allAtoms.bonds[0]))
        

    def set_amide_torsions(self, flag):
        if not len(self.molecule.amidebnds):
            return
        self.RBM.set_amide_torsions(flag)
        
        
    def set_ppbb_torsions(self, flag):
        if not len(self.molecule.ppbbbnds):
            return
        self.RBM.set_peptidebackbone_torsions(flag)
        
        
    def set_guanidinium_torsions(self, flag):
        if not len(self.molecule.guanidiniumbnds):
            return
        self.RBM.set_guanidinium_torsions(flag)
        

    def toggle_torsion(self,ind1, ind2):
        self.RBM.toggle_torsion(ind1, ind2)
        

    def changePlanarityCriteria(self, cutoff):
        oldAromCs = self.aromCs
        self.aromCs = self.ACM.setAromaticCarbons(self.molecule, cutoff)
        if self.debug: print "now there are ", len(self.aromCs), " aromaticCs"


    def set_carbon_names(self, atoms, type):
        self.ACM.set_carbon_names(atoms, type)
        #have to reset this field here
        self.aromCs = AtomSet(filter(lambda x: (x.autodock_element=='A' and x.element=='C'),
                                          self.molecule.allAtoms))


    def setAromaticCarbons(self):
        self.aromCs = self.ACM.setAromaticCarbons(self.molecule)
        return self.aromCs

#LigandPreparation

class AD4LigandPreparation(LigandPreparation):

    def __init__(self, molecule, mode='automatic', repairs='checkhydrogens',
                    charges_to_add='gasteiger', cleanup='nphs_lps_waters_nonstdres',
                    allowed_bonds='backbone',  #backbone on, amide + guanidinium off
                    root='auto', outputfilename=None, dict=None, debug=False,
                    check_for_fragments=False, bonds_to_inactivate=[],
                    inactivate_all_torsions=False,
                    version=4, typeAtoms=True, limit_torsion_number=False, 
                    limit_torsions_type=None, limit_torsions=False,
                    delete_single_nonstd_residues=False,
                    attach_nonbonded_fragments=False,
                    detect_bonds_between_cycles=False):



        if attach_nonbonded_fragments==True:
            molecule.attach_nonbonded_fragments()
        #FIX THIS: what if molecule already has autodock_element set???
        LigandPreparation.__init__(self, molecule, mode='interactive', 
                    repairs=repairs, charges_to_add=charges_to_add, 
                    cleanup=cleanup, allowed_bonds=allowed_bonds, 
                    root=root, outputfilename=outputfilename, 
                    dict=dict, debug=debug,
                    check_for_fragments=check_for_fragments,
                    bonds_to_inactivate=bonds_to_inactivate, 
                    inactivate_all_torsions=inactivate_all_torsions,
                    version=version,
                    limit_torsions=limit_torsions,
                    detect_bonds_between_cycles=detect_bonds_between_cycles)

        # AD4 force field uses total rotatable bonds for TORSDOF:
        # calculation of loss of entropy on ligand binding
        #molecule.TORSDOF = molecule.possible_tors
        self.charge_error_tolerance = 0.008
        if typeAtoms:
            #NEED TO TYPE ATOMS!!
            delattr(molecule.allAtoms, 'autodock_element')
            atomTyper = AutoDock4_AtomTyper()
            atomTyper.setAutoDockElements(molecule, reassign=True) #4/15/05:catch here?
        #if self.debug: print "allAtoms.autodock_element=", molecule.allAtoms.autodock_element
        molecule.TORSDOF = molecule.possible_tors - len(molecule.guanidiniumbnds+molecule.amidebnds)
        #could limit the number of torsions here
        if limit_torsion_number is not False:
            ndihe = limit_torsion_number
            type = limit_torsions_type
            self.RBM.limit_torsions(ndhihe, type)

        self.writer = AD4LigandWriter()
        #MODE SWITCH 5: write outputfile     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: what records should be written?
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def set_autodock_element_types_for_writing(self, allAtoms):
        pass 

    def restore_autodock_element_types_after_writing(self, problem_atoms):
        pass



class LigandWriter:

    """
    LigandWriter class writes a ligand molecule which has a LigandPreparationObject to a file
    """

    def __init__(self, debug=False, write_CONECT=False):
        self.writer = PdbqWriter()
        self.maxtors = 32
        self.debug = debug
        self.write_CONECT = write_CONECT


    def write(self, ligand, outputfile):
        #setup rnum field
        ligand.allAtoms.rnum = -1

        old_torscount=ligand.torscount
        ligand.torscount = len(ligand.allAtoms.bonds[0].get(lambda x:x.activeTors))
        if self.debug and old_torscount!=ligand.torscount:
            print "reset torscount from ", old_torscount, " to ", ligand.torscount
        
        ttc = ligand.torscount
            
        assert hasattr(ligand, 'ROOT') # is this necessary or redundant
        assert hasattr(ligand, 'TORSDOF')
        self.outfptr = open(outputfile, 'w')
        #FIRST write out remarks about torsions
        if (ttc>self.maxtors):
            self.outfptr.write("REMARK WARNING: %d MAX_TORS EXCEEDED!!!\n"%(self.maxtors))
        self.outfptr.write("REMARK  " +"%d" %(ttc) + " active torsions:\n")
        self.outfptr.write("REMARK  status: ('A' for Active; 'I' for Inactive)\n")
        self.outatom_counter = 1
        bondset = ligand.allAtoms.bonds[0]
        rootlist = [ligand.ROOT]
        rootnum = 1
        for b in bondset:
            if b.activeTors==1: 
                bTors = 'A'
            else:   
                bTors = 'I'
            at1str = b.atom1.name + '_' + str(b.atom1.number)
            at2str = b.atom2.name + '_' + str(b.atom2.number)
            if b.possibleTors :
                if bTors=='A' :
                    outstring = "REMARK " +" %3d  %s    between atoms: %-3s  and  %-3s \n" %(self.outatom_counter,bTors, at1str, at2str)
                    b.outNumber = self.outatom_counter
                    self.outatom_counter = self.outatom_counter + 1
                else:
                    outstring = "REMARK " +"      %s    between atoms: %-2s  and  %-3s \n" %(bTors, at1str, at2str)
                    #outstring = "REMARK " +"      %s    between atoms: %-3s  and  %-3s \n" %(bTors,b.atom1.name,b.atom2.name)
                self.outfptr.write(outstring)
        #next write out  root
        self.outfptr.write("ROOT\n")
        self.outatom_counter = 1
        #reset used field to serve as visited flag
        ligand.allAtoms.used = 0
        #rootlist grows to include atoms up to first active tors in each subtree
        #this starts w/ 1 atom:
        self.writtenAtoms = []
        for at in rootlist:
            #to start, rootlist has just 1 atom
            at.used = 1
            for bond in at.bonds:
                #if bond.activeTors and bond.possibleTors:continue
                at2 = bond.atom1
                if at2==at: 
                    at2 = bond.atom2
                if bond.activeTors and bond.possibleTors:
                    #if at2 not in branchAtoms:
                        #branchAtoms.append(at2)
                    continue
                if at2.used: 
                    continue
                if at2 not in rootlist:
                    rootlist.append(at2)
                    if self.debug: print 'added ', at2.name, ' to rootlist'
                    at2.rnum0 = rootnum
                    rootnum = rootnum + 1
        #before writing, sort rootlist
        newrootlist = []
        oldrootlist = rootlist
        r_ctr = 0
        for at in ligand.chains.residues.atoms:
            if hasattr(at, 'rnum0'):
                newrootlist.append(at)
                at.rnum = r_ctr
                r_ctr = r_ctr+1
        for at in newrootlist:
            at.number = self.outatom_counter
            self.writer.write_atom(self.outfptr,at)
            at.newindex = self.outatom_counter
            at.used = 1
            self.writtenAtoms.append(at)
            self.outatom_counter = self.outatom_counter+1
        #at this point, add chain_roots if there are any
        #for ch in dict['chain_rootlist']:
        #    for at in ch.residues.atoms:
        #        at.number = self.outatom_counter
        #        at.rnum = r_ctr
        #        r_ctr = r_ctr+1
        #        self.writer.write_atom(self.outfptr,at)
        #        at.newindex = self.outatom_counter
        #        at.used = 1
        #        self.writtenAtoms.append(at)
        #        self.outatom_counter = self.outatom_counter+1
        self.outfptr.write("ENDROOT\n")
        #last write out the rest of the stuff, using WriteSubtree.....
        #for at in rootlist:
        for at in rootlist:
            for bond in at.bonds:
                at2 = bond.atom1
                if at2==at: 
                    at2 = bond.atom2
                if at2.used: 
                    continue
                self.process(at, at2)
###                marker = self.outatom_counter
###                outstring = "BRANCH %3d %3d\n"%(at.rnum+1, marker)
###                self.outfptr.write(outstring)
###                ###11/30: change #1-> write atom before calling WriteSubtree
###                at2.newindex = self.outatom_counter
###                at2.number = self.outatom_counter
###                self.writer.write_atom(self.outfptr,at2)
###                self.writtenAtoms.append(at2)
###                self.outatom_counter = self.outatom_counter+1
###                self.WriteSubtree(at,at2)
###                outstring = "ENDBRANCH %3d %3d\n"%(at.rnum +1,marker)
###                self.outfptr.write(outstring)
        outstring = "TORSDOF " + str(ligand.TORSDOF) + "\n"
        self.outfptr.write(outstring)
        #call write_CONECT here
        if self.write_CONECT:
            self.writeCONECTRecords(ligand, self.outfptr)
        self.outfptr.close()
        ligand.returnCode = 0

        if len(self.writtenAtoms)!=len(ligand.allAtoms):
            ligand.returnCode = 1
            allAtoms = len(ligand.allAtoms)
            notWritten = allAtoms - len(self.writtenAtoms)
            ligand.returnMsg = 'WARNING: %d atoms of %d in %s  were not written\n' %(notWritten, allAtoms, ligand.parser.filename)
        for a in ligand.allAtoms:
            if hasattr(a, 'rnum0'):
                delattr(a, 'rnum0')
            if hasattr(a, 'cycleout'):
                delattr(a, 'cycleout')
            if hasattr(a, 'newindex'):
                a._ni = a.newindex
                delattr(a, 'newindex')
        ligand.ROOT.rnum0 = 0
        ligand.outputfile = outputfile


    def process(self, fromAtom, nextAtom):
        startIndex = fromAtom.number
        endIndex = self.outatom_counter
        outstring = "BRANCH %3d %3d\n"%(startIndex, endIndex)
        #outstring = "BRANCH %3d %3d\n"%(fromAtom.rnum+1, marker)
        self.outfptr.write(outstring)
        queue = self.writeBreadthFirst(self.outfptr, fromAtom, nextAtom)
        if self.debug: print fromAtom.name, ':', nextAtom.name, ': queue=', queue
        if len(queue):
            for fromAtom, nextAtom in queue:
                if self.debug: print " processing queue entry: ", fromAtom.name, '-', nextAtom.name
                self.process(fromAtom, nextAtom)
        outstring = "ENDBRANCH %3d %3d\n"%(startIndex, endIndex)
        self.outfptr.write(outstring)


    def writeLevel(self, atom, outfptr):
        """
        write all atoms bonded to atoms bonded to this atom by non-rotatable
        bonds
        """
        if self.debug: 
            print "\n\nin writeLevel with ", atom.name, " outatom_counter=", self.outatom_counter
            print "len(", atom.name, ").bonds=", len(atom.bonds)
        queue = []
        nextAts = []
        for b in atom.bonds:
            if self.debug:
                print "processing b=", b.atom1.name, '-', b.atom2.name, ' activeTors=', b.activeTors
                print 'atom1 in writtenAtoms=', b.atom1 in self.writtenAtoms
                print 'atom2 in writtenAtoms=', b.atom2 in self.writtenAtoms
            if b.activeTors: 
                at2 = b.atom1
                if at2==atom: at2=b.atom2
                queue.append((atom, at2))
                if self.debug: print atom.name, 'wL: queue=', queue
                continue
            a2 = b.atom1
            if a2==atom:
                a2 = b.atom2
            if a2.used: 
                if self.debug: print "!!a2 is already used!!", a2.name
                continue    
            if a2 not in self.writtenAtoms:
                a2.number = a2.newindex = self.outatom_counter
                if self.debug: print "writeLevel: wrote bonded atom named=", a2.name, 'a2.used=', a2.used
                self.writer.write_atom(outfptr, a2)
                self.writtenAtoms.append(a2)
                a2.used = 1
                self.outatom_counter+=1
                nextAts.append(a2)
        for a2 in nextAts:
            if self.debug: 
                print 'in for nextAts loop with a2=', a2.name
                print 'calling wL'
            nq = self.writeLevel(a2, outfptr)
            if len(nq):
                if self.debug: print "extending queue with", nq
                queue.extend(nq)
        if self.debug:
            print 'returning queue=', queue
        return queue
            

    def writeBreadthFirst(self, outfptr, fromAtom, startAtom):
        """
            None <-writeBreadthFirst(outfptr, fromAtom, startAtom)
            writeBreadthFirst visits all the atoms in the current level
            then the first level down etc in a Breadth First Order traversal. 
                            1                <-1
                        5 6   7 8            <-3
                     9 10   11 12            <-4
            It is used to write out the molecule with the correct format 
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH 
            statements are added. 
        """
        if self.debug: 
            print "in wBF with fromAtom=", fromAtom.name, '+ startAtom=', startAtom.name, 'startAtom.used=', startAtom.used
        queue = []
        if startAtom.used==0:
            startAtom.used = 1
            startAtom.number = startAtom.newindex = self.outatom_counter
            self.writer.write_atom(outfptr,startAtom)
            if self.debug: print 'wBF: wrote ', startAtom.name
            self.writtenAtoms.append(startAtom)
            self.outatom_counter += 1
            if self.debug: print "self.outatom_counter=", self.outatom_counter
            activeTors = []
            #outfptr.write(outstring)
            for bond in startAtom.bonds:
                at2 = bond.atom1
                if at2==startAtom: at2 = bond.atom2
                if at2==fromAtom: continue  #this is current bond
                elif not at2.used:
                    if bond.activeTors:
                        queue.append((startAtom,at2))
                    else:
                        at2.number = at2.newindex = self.outatom_counter
                        if self.debug: 
                            print "\n\nwriting and calling wL with nA=", at2.name, '-', at2.number
                        self.writer.write_atom(outfptr, at2)
                        if self.debug: print 'wBF2: wrote ', at2.name
                        at2.written = 1
                        self.writtenAtoms.append(at2)
                        at2.newindex = self.outatom_counter
                        self.outatom_counter = self.outatom_counter + 1
                        if self.debug: print '!!!2:calling wL'
                        newQ = self.writeLevel(at2, outfptr)
                        if self.debug: print "newQ=", newQ
                        at2.used = 1
                        if len(newQ): 
                            if self.debug: print "@@@@len(newq)=", len(newQ)
                            queue.extend(newQ)
                            if self.debug: print "queue=", queue
            if self.debug: 
                print " currently queue=",
                for atom1, atom2 in queue: 
                    print atom1.name, '-', atom2.name, ',',
                print
        return  queue


    def WriteSubtree(self,fromAtom, startAtom):
        """WriteSubtree recursively visits the atoms in the current 'subtree' of the molecule
        in a BREADTH First Order traversal. It is used to write out the molecule with
        the correct format for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH and
        TORS/ENDTORS statements are added. 
        """
        
        if startAtom.used==0:
            startAtom.used = 1
            at = startAtom
            for bond in startAtom.bonds:
                if bond.activeTors: 
                    continue
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom: 
                    nextAtom = bond.atom2
                if nextAtom==fromAtom: 
                    continue
                if not nextAtom.used:
                    if hasattr(bond,'incycle'):
                        if not hasattr(nextAtom, 'cycleout'):
                            nextAtom.cycleout = 1
                            nextAtom.newindex = self.outatom_counter
                            nextAtom.number = self.outatom_counter
                            self.writer.write_atom(self.outfptr,nextAtom)
                            self.writtenAtoms.append(nextAtom)
                            self.outatom_counter = self.outatom_counter+1
                    else:
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(self.outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
            for bond in startAtom.bonds:
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom: 
                    nextAtom = bond.atom2
                if nextAtom==fromAtom: 
                    continue
                if not nextAtom.used:
                    testcond = len(nextAtom.bonds)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "BRANCH %3d %3d\n"%(at.newindex,marker)
                            self.outfptr.write(outstring)
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(self.outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
                    self.WriteSubtree(startAtom, nextAtom)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "ENDBRANCH %3d %3d\n"%(at.newindex,marker)
                            self.outfptr.write(outstring)
        return


    def writeCONECTRecords(self, mol, fptr):
        atms = mol.allAtoms
        for atm in atms:
            rec = 'CONECT%5i'%atm.number
            for b in atm.bonds:
                a2 = b.atom1
                if a2==atm: a2 = b.atom2
                if not a2 in atms: continue #don't write intermolecular bonds
                rec = rec + '%5i'%a2.number 
            rec = rec + '\n'
            fptr.write(rec)
            


class AD4LigandWriter(LigandWriter):
    """
    AD4LigandWriter class writes a ligand molecule which has a LigandPreparationObject to a file
    """

    def __init__(self, debug=False, write_CONECT=False):
        LigandWriter.__init__(self, write_CONECT)
        self.writer = PdbqtWriter()
        self.maxtors = 32
        self.debug = debug




class RotatableBondManager:
    """
    flags bonds of ligand molecules for possible rotability and active 
    rotability.  [for backwards compatibility] The flags are 'possibleTors' 
    and 'activeTors'
    """


    def __init__(self, molecule, allowed_bonds,
                    root='auto',  debug=False, 
                    check_for_fragments=False,
                    bonds_to_inactivate='',
                    detectAll=False):
        if debug: 
            print "bonds_to_inactivate=", bonds_to_inactivate, " with len=", 
            print len(bonds_to_inactivate)
        self.detectAll = detectAll
        allowed_bond_list = []
        if allowed_bonds:
            allowed_bond_list = string.split(allowed_bonds,'_')
        #set up flags:
        allow_amide_torsions = 'amide' in allowed_bond_list
        molecule.has_amide = allow_amide_torsions
        if debug: print "allow amide=", allow_amide_torsions
        allow_backbone_torsions = 'backbone' in allowed_bond_list
        molecule.has_backbone = allow_backbone_torsions
        allow_guanidinium_torsions = 'guanidinium' in allowed_bond_list
        molecule.has_guanidinium = allow_guanidinium_torsions
        allow_all_torsions = 'all' in allowed_bond_list
        self.molecule = molecule
        self.debug = debug
        self.__classifyBonds(molecule.allAtoms, allow_guanidinium_torsions)
        self.set_peptidebackbone_torsions(allow_backbone_torsions)
        self.set_amide_torsions(allow_amide_torsions)
        #self.set_guanidinium_torsions(allow_guanidinium_torsions)
        #NEW 11/22/2004
        if check_for_fragments:
            self.detect_bonded_fragments()
            
        if root=='auto':
            self.autoroot()
        else:
            self.setroot(int(root))
        if len(bonds_to_inactivate):
            bnds_list = map(int,string.split(bonds_to_inactivate,'_'))
            #???
            #molecule.has_amide = 'amide' not in bnds_list
            #molecule.has_guanidinium = 'guanidinium' not in bnds_list
            #molecule.has_backbone = 'backbone' not in bnds_list
            #molecule.has_active = 'all' not in bnds_list
            if debug: print "bnds_list=", bnds_list
            for i in range(len(bnds_list)/2):
                ind1 = bnds_list[i*2]
                ind2 = bnds_list[i*2+1]
                if debug: print "inactivating ", ind1, ' to ', ind2
                self.toggle_torsion(ind1, ind2)


    def __classifyBonds(self, atoms, allow_guanidinium_torsions):
        mols = atoms.top.uniq()
        #check that all are in the same molecule, for now
        assert  len(mols)==1
        mol = mols[0]
        if hasattr(mol,'processed_bonds'):
            return
        #check whether these atoms already have bonds
        if not len(atoms.bonds[0]):
            mol.buildBondsByDistance()
        ADBC = self.ADBC = AutoDockBondClassifier(detectAll=self.detectAll)
        dict =self.dict = ADBC.classify(mol.allAtoms.bonds[0])
        #turn everything off 
        mol.allAtoms.bonds[0].possibleTors = 0
        mol.allAtoms.bonds[0].activeTors = 0
        mol.allAtoms.bonds[0].leaf = 0
        mol.allAtoms.bonds[0].incycle = 0
        # restore appropriate categories:
        if len(dict['leaf']): 
            dict['leaf'].leaf = 1
        aromBnds = dict['aromatic']
        if len(dict['cycle']): 
            dict['cycle'].incycle = 1
        if len(dict['rotatable']): 
            if self.debug: print "len(dict['rotatable'])=", len(dict['rotatable'])
            for b in dict['rotatable']:
                #check for mistake
                if b.incycle: 
                    print "removing rotatable ERROR: bnd", b
                    b.possibleTors = 0
                    b.activeTors = 0
                    dict['rotatable'].remove(b)
            dict['rotatable'].possibleTors = 1
            dict['rotatable'].activeTors = 1
            #NB: this should set amide and ppbb
        #mol.guanidiniumbnds = BondSet()
        #self.set_guanidinium_torsions(allow_guanidinium_torsions)
        for b in dict['guanidinium']:
            if self.debug: print b, ' ', b.bondOrder
            if not b.incycle and not b.leaf and b.bondOrder==1:
                if self.debug: print 'set ', b,'.possible/active to 1/0'
                b.possibleTors = 1 
                b.activeTors = allow_guanidinium_torsions
        #if len(dict['guanidinium']): and not allow_guanidinium_bonds: 
        #    #as per gmmm 4/14:
        #    mol.guanidiniumbnds = dict['guanidinium']
        #    for b in dict['guanidinium']:
        #        if self.debug: print b, ' ', b.bondOrder
        #        if not b.incycle and not b.leaf and b.bondOrder==1:
        #            if self.debug: print 'set ', b,'.possible/active to 1/0'
        #            b.possibleTors = 1 
        #            b.activeTors = 0
        #set keywords for this molecule
        #NB: 'torscount' must be adjusted if specific torsions are inactivated
        rotatable_bnds = mol.allAtoms.bonds[0].get(lambda x: x.activeTors==1)
        possible_tors_bnds = mol.allAtoms.bonds[0].get(lambda x: x.possibleTors==1)
        possible_tors_ct = 0
        mol.possible_tors_bnds = BondSet()
        if possible_tors_bnds:
            possible_tors_ct = len(possible_tors_bnds)
            mol.possible_tors_bnds = possible_tors_bnds
        mol.torscount = 0
        mol.possible_tors = possible_tors_ct
        if rotatable_bnds:
            #THIS IS WRONG!!rotatable includes too many
            #mol.torscount = len(rotatable_bnds)
            #torscount changes, possible_tors_ct doesn't
            mol.torscount = possible_tors_ct
        if self.debug:
            for b in rotatable_bnds:
                if b not in dict['rotatable']:
                    print "ERROR in rotatable_bnds set:"
                    print b.atom1.full_name(), '-', b.atom2.full_name()
        if self.debug: print "set ", mol.name, ".torscount=", mol.torscount
        #mol.TORSDOF = mol.torscount - len(dict['hydrogenRotators'])
        mol.TORSDOF = possible_tors_ct - len(dict['hydrogenRotators'])
        mol.possible_tors = possible_tors_ct
        if self.debug:
            for b in dict['hydrogenRotators']:
                print "possible_tors_ct =", possible_tors_ct
                print "len(dict['hydrogenRotators']) =", len(dict['hydrogenRotators'])
                print "set ", mol.name, ".TORSDOF=", mol.TORSDOF
        mol.amidebnds = dict['amide']
        mol.ppbbbnds = dict['ppbb']
        mol.guanidiniumbnds = dict['guanidinium']
        mol.processed_bonds = 1


    def set_torsions(self, mol, bnds, flag):
        ct = 0
        for b in bnds:
            if not b.possibleTors:
                continue
            if flag:
                if not b.activeTors:
                    b.activeTors = 1
                    mol.torscount = mol.torscount + 1
                    ct = ct + 1
            else:
                if b.activeTors==1:
                    b.activeTors = 0
                    mol.torscount = mol.torscount - 1
                    ct = ct - 1
        return ct


    def set_amide_torsions(self, flag): 
        """ 
        set rotability of all amide bonds to flag: True/False
        """
        if self.debug: print "in set_amide_torsions with flag", flag
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        #print "before call to set_torsions: mol.amidebnds.activeTors=", mol.amidebnds.activeTors
        ct = self.set_torsions(mol, mol.amidebnds, flag)
        if self.debug: 
            print mol.name, ':', ct, " amide torsions: ", mol.name, 
            print ".torscount=", mol.torscount


    def set_peptidebackbone_torsions(self, flag): 
        """ 
        set rotability of all peptidebackbone bonds to flag: True/False
        """
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        ct = self.set_torsions(mol, mol.ppbbbnds, flag)
        if self.debug: 
            print mol.name,':',ct, " peptidebackbone torsions: ", mol.name,
            print ".torscount=", mol.torscount


    def set_guanidinium_torsions(self, flag): 
        """ 
        set rotability of all guanidinium bonds to flag: True/False
        """
        if self.debug: print "SGT: flag=", flag
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        ct = self.set_torsions(mol, mol.guanidiniumbnds, flag)
        if self.debug: 
            print mol.name,':',ct, " guanidinium torsions: ", mol.name,
            print ".torscount=", mol.torscount


    def set_all_torsions(self, flag): 
        """ 
        set rotability of all rotatable bonds to flag: True/False
        """
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        rotatable = filter(lambda x: x.possibleTors==1, mol.allAtoms.bonds[0])
        if len(rotatable):
            ct = self.set_torsions(mol, rotatable, flag)
            if self.debug: 
                print mol.name, ':', ct, " active torsions: ", mol.name, 
                print ".torscount=", mol.torscount


    def toggle_torsion(self, ind1, ind2, use_tt_ind=0):
        """ 
        flip the rotability of the bond 
        between atoms with indices ind1 and ind2
        use_tt_ind is set when limiting torsion using a built-on-the-fly
        rotabilityManager
        """
        mol = self.molecule
        at1 = mol.allAtoms[ind1]
        at2 = mol.allAtoms[ind2]
        if use_tt_ind:
            at1 = mol.allAtoms.get(lambda x: x.tt_ind==ind1)
            if not at1: 
                print "unable to toggle torsion with tt_ind[0]=", ind1
                return "ERROR"
            at1 = at1[0]
            at2 = mol.allAtoms.get(lambda x: x.tt_ind==ind2)
            if not at2: 
                print "unable to toggle torsion with tt_ind[1]=", ind2
                return "ERROR"
            at2 = at2[0]
        bnds = AtomSet([at1,at2]).bonds[0]
        if not len(bnds):
            print "ERROR: atoms at indices ", ind1, ' and ', ind2, ' have no common bond'
            return 
        bnd = bnds[0]
        if bnd.possibleTors:
            #print "calling set_torsions with ", bnd, " and flag=", abs(bnd.activeTors-1)
            self.set_torsions(mol, [bnd], abs(bnd.activeTors-1))
        else:
            print "ERROR: bond between these atoms is non-rotatable"
            return 

        
    def detect_bonded_fragments(self):
        atoms = self.molecule.allAtoms
        #first: build a list of bonded fragments:
        l = []   # the list to hold dicts of bonded fragments
        fb = []  # keep track of which bonds are found
        #start with unique list of bonds 'bnds'
        bond_dict = {}
        for b in atoms.bonds[0]:
            bond_dict[b] = 1  
        bnds = bond_dict.keys()
        #first pass: try to connect pieces
        for b in bnds:
            ind1 = b.atom1
            ind2 = b.atom2
            found = b in fb
            if not found:
                for d in l:
                    d_keys = d.keys()
                    if ind1 in d_keys:
                        d[ind2] = 1
                        fb.append(b)
                    elif ind2 in d_keys:
                        d[ind1] = 1
                        fb.append(b)
                if not b in fb:
                    #start new fragment 
                    l.append({ind1:1, ind2:1})
                    fb.append(b)
        #now to try to merge fragments, use lists of keys
        key_list = []
        for dict in l:
            key_list.append(dict.keys())
        #check each list of keys against following lists
        #if there are any duplications, merge current
        #into the following one..
        for i in range(len(key_list)):
            #check this list 
            kl = key_list[i]
            found = 0
            #...against the each of the subsequent ones
            for j in range(i+1, len(key_list)):
                jl = key_list[j]
                #....check each entry in subsequent one
                #.........against this list 'kl'
                for entry in jl:
                    if entry in kl:
                        #if a match
                        #...merge this dict into jth one..
                        l[j].update(l[i])
                        #.....reset this dict to {}
                        l[i] = {}
                        #.......set found flag
                        found = 1
                        #..........update jth list of keys
                        key_list[j]= l[j].keys()
                        #............skip rest of jl
                        break
                if found:
                    #.................and skip rest of key_list
                    break
        #now find an index in the largest fragment
        max_ind = 0
        largest_fragment = {}
        ct = 0
        for dict in l:
            len_d = len(dict)
            if len_d >0:
                ct = ct + 1
            if len_d > max_ind:
                max_ind = len(dict)
                largest_fragment = dict
                #the keys are ATOMS!
                #print z_keys
        if ct>0:
            #build an atomset of the keys
            self.molecule.largest_fragment = AtomSet(largest_fragment.keys())
        
    
    def autoroot(self, verbose=0):
        """
        autoRoot
        """
        mol = self.molecule
        #clear old root
        if hasattr(mol, 'ROOT') and hasattr(mol.ROOT, 'rnum0'):
            if verbose: print "delattr mol.ROOT.rnum0"
            delattr(mol.ROOT, 'rnum0')
        if hasattr(mol, 'autoRoot'):
            mol.ROOT = mol.autoRoot
            mol.ROOT.rnum0 = 0
            if verbose: print "mol.ROOT exists! setting ROOT.rnum0"
            return
        if len(mol.chains)>1:
            return  "AutoRoot not implemented for molecules with >1 chain"
        if hasattr(mol, 'largest_fragment'):
            if verbose: print "mol.largest_fragment exists=", mol.largest_fragment 
            atoms_to_check = mol.largest_fragment
        else:
            atoms_to_check = mol.allAtoms
            if verbose: print "mol.largest_fragment doesn't exists=> using ", len(mol.allAtoms)
        #mol.bestBranch = len(mol.allAtoms)
        mol.bestBranch = atoms_to_check
        if verbose: print "length mol.bestBranch=", len(mol.bestBranch)
        bestList = []
        #for item in mol.allAtoms:
        for item in atoms_to_check:
            if not hasattr(item, 'bonds'): 
                if verbose: print 'skip ', item.number, ': no bonds'
                continue
            if len(item.bonds)==1 and item.bonds[0].leaf: 
                if verbose: print 'skip ', item.number, ': leaf atom'
                continue
            if hasattr(item,'leaf'): 
                if verbose: print 'skip ', item.number, ': leaf atom'
                continue
            if verbose: print "processing atom ", item.number
            item.maxbranch = 0
            for b in item.bonds:
                nxtAtom = b.atom1
                if nxtAtom==item:
                    nxtAtom = b.atom2
                if not b.leaf:
                    thistree = mol.subTree(item, nxtAtom, atoms_to_check)
                    #thistree = mol.subTree(item, nxtAtom, mol.allAtoms)
                    thisbranch = len(thistree)
                    if verbose: print 'subtree for bond ', b.atom1.number, '-', b.atom2.number, ' consists of ', thisbranch, ' atoms'
                    if thisbranch>item.maxbranch:
                        if verbose: print item.number, '-', item.name,': new largest subtree=', thisbranch
                        item.maxbranch = thisbranch
            #bestList holds all current best choices for Root..
            if item.maxbranch<mol.bestBranch:
                if verbose: print 'NEW ROOT CANDIDATE ', item.number,'-', item.name, ': largest subtree length =', item.maxbranch
                bestList = []
                bestList.append(item)
                mol.bestBranch = item.maxbranch
                if verbose: 
                    print 'currently there is/are ', len(bestList), ' root candidates:'
                    for a in bestList: 
                        print a.full_name(),
                    print
            if item.maxbranch==mol.bestBranch and item not in bestList:
                bestList.append(item)
                if verbose: print 'added item', len(bestList)
        if len(bestList)>1:
            foundCycle = 0
            for at in bestList:
                at.cycleatom = 0
                for b in at.bonds:
                    if hasattr(b, 'incycle') and b.incycle:
                        at.cycleatom = 1
                        continue
            for at in bestList:
                if at.cycleatom:
                    mol.ROOT = at
                    mol.autoRoot = at
                    mol.ROOT.rnum0 =0
                    foundCycle = 1
                    break
            #if bestList had a cycle atom, it's been set to root..if NOT:
            if verbose: print 'foundCycle= ', foundCycle
            if not foundCycle:
                mol.autoRoot = bestList[0]
                mol.ROOT = bestList[0]
                mol.ROOT.rnum0 =0
        #if ties for possible root, use first entry in bestRoot list...
        elif len(bestList):
            if verbose: 
                print 'using first entry in bestList of ', len(bestList)
                for item in bestList:
                    print item.name,
                print
            mol.autoRoot = bestList[0]
            mol.ROOT = bestList[0]
            mol.ROOT.rnum0 =0
        else:
            #if the list is empty, use the first atom
            # this is a correction added for 'HF' ligand
            if verbose: print 'no entries in bestList using first atom '
            mol.autoRoot = mol.allAtoms[0]
            mol.ROOT = mol.autoRoot
            mol.ROOT.rnum0 =0
        return mol.ROOT


    def setroot(self, index):
        """
        setroot to atom at index
        """
        #if autoRoot has already been set, remove its 'rnum0'
        mol = self.molecule
        if hasattr(mol, 'autoRoot') and hasattr(mol.autoRoot, 'rnum0'):
            delattr(mol.autoRoot, 'rnum0')
        if hasattr(mol, 'ROOT') and hasattr(mol.ROOT, 'rnum0'):
            if self.debug:
                print "deleting ", mol.ROOT.name, ' rnum0'
            delattr(mol.ROOT, 'rnum0')
        #set the root to indicated atom
        mol.ROOT = mol.allAtoms[index]
        #set the root's rnum0 field
        mol.ROOT.rnum0 = 0
        return


    def limit_torsions(self, numTors, type='fewest'):
        """
        numTors, type='fewest'
        limit torsions to a specified number, numTors
        #changed per GMM request: invert logic of type here
        previously 'type' meant TOGGLE torsions moving 'type' atoms: 
        now 'type' means KEEP torsions moving type atoms: 
        options are 
                'fewest' atoms which is the default
              or 'most' 
        """
        mol = self.molecule
        if not hasattr(mol, 'ROOT'):
            print 'no ROOT specified for ', mol.name, 
            print ' unable to limit torsions'
            return 'ERROR'
        if not hasattr(mol, 'torTree') or \
                not hasattr(mol.allAtoms[0], 'tt_ind'):
            mol.torTree = TorTree(mol.parser, mol.ROOT)
        self._setTorsions(numTors, type)


    def _setTorsions(self, numTors, type):
        #print "in _setTorsions ", numTors, " and type ", type
        assert type in ['fewest', 'most'], 'unknown torsion type'
        mol = self.molecule
        torsionMap = mol.torTree.torsionMap
        tNum = len(torsionMap)
        if numTors>tNum:
            print 'too many torsions specified! '+ str(numTors)+  ' reducing to'+str(tNum)
            numTors = tNum
        #set up rangeList to implement which kind of torsion
        # to turn off or on
        #fewest uses range 0,1,2,.... (from beginning toward end)
        #most uses range -1, -2,-3,.... (from end toward beginning)
        if type=='fewest':
            rangeList = range(numTors)
        else:
            rangeList = []
            for k in range(1, numTors+1):
                rangeList.append(-k)
        #turn them all off
        ct = self.set_torsions(mol, mol.allAtoms.bonds[0], 0)
        #print "ct = ", ct, " and mol.torscount=", mol.torscount
        #torsionMap = mol.torTree.torsionMap
        #for i in range(tNum):
        #    node = torsionMap[i]
        #    ind1, ind2 = node.bond
        #    self.toggle_torsion(ind1, ind2, use_tt_ind=1)

            #b = mol.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            #b.activeTors = 0
            #mol.torscount -= 1
        #turn on the right number at correct end
        for i in rangeList:
            node = torsionMap[i]
            ind1, ind2 = node.bond
            self.toggle_torsion(ind1, ind2, use_tt_ind=1)
            #b = mol.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            #b.activeTors = 1
            #mol.torscount += 1



# DEPRECATED 2009
class AD4FlexibleDockingPreparation:
    """
    DEPRECATED 2009
    Facade for preparing input for an AutoDock4 experiment 
    with ligand plus flexible sidechains in the receptor in one file

    This preparation involves:
        1. adding autodock_element types 
        2. selecting flexible residues 
        3. specifying patterns of flexibility in side-chains of selected
        residues
        4. specifying ligand 
        5. writing a flexible pdbqt outputfile containing ligand plus
        side-chains
        6. writing a rigid pdbqt outputfile containing remaining receptor
        atoms

    optionally
        writing files for a directory of ligands with a single receptor

    Receptor preparation collaborators extend baseclass collaborators:
        in AutoDockTools:
            in atomTypeManager classes for conforming to AD atomtypes:
                AutoDock4_AtomTyper
       in this file:
            AD4ReceptorWriter
        
"""

    def __init__(self, molecule, ligand, mode='automatic',
                    rigid_filename=None,  residues=[],
                    allowed_bonds='',  #backbone on, amide + guanidinium off
                    non_rotatable_bonds=[], flex_filename=None,
                    reassign=True, debug=False):



        #?NEED TO TYPE ATOMS?
        self.molecule = molecule
        self.debug = debug
        file_type = os.path.splitext(os.path.basename(molecule.parser.filename))[1]
        if file_type=='.pdbqt':
            pass
        elif file_type in ['.pdbqs','.pdbq','.mol2', '.pdb']:
            if file_type in ['.pdbqs','.pdbq']:
                try:
                    delattr(molecule.allAtoms, 'autodock_element')
                except:
                    pass
            atomTyper = AutoDock4_AtomTyper()
            atomTyper.setAutoDockElements(molecule, reassign=reassign) #
        else:
            print file_type, " unknown type: must be .pdbqt, .pdbqs, .pdbq, .pdb or .mol2" 
        self.ligand = ligand
        self.writer = PdbqtWriter()
        if not hasattr(ligand, 'torTree'):
            if self.debug: print "creating ligand.lpo"
            #???which other input options should be supported???
            self.ligand.lpo = AD4LigandPreparation(ligand, charges_to_add='gasteiger')
        try:
            self.allowed_bonds = allowed_bonds.split('_')  # a string "backbone"
        except:
            self.allowed_bonds = ['']
        self.non_rotatable_bonds = non_rotatable_bonds
        # have to get the number of torsions in the residues before writing
        # the ligand to the outputfilename
        flex_residues = self.flex_residues = self.process_residues(residues)
        #MODE SWITCH 5: write outputfiles     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: are there other records which should be written?
            if flex_filename is None:
                flex_filename = self.ligand.name+'_'+self.molecule.name+'_flex.pdbqt'
            if rigid_filename is None:
                rigid_filename = self.ligand.name+'_'+self.molecule.name+'_rigid.pdbqt'
            self.write_flex(flex_residues, self.ligand.parser.filename, flex_filename)
            if self.debug: print "wrote flexible file:", flex_filename
            self.write_rigid(self.molecule, rigid_filename)
            if self.debug: print "wrote rigid file:", rigid_filename


    def write(self, outputfilenames=None):
        if outputfilenames is not None and len(outputfilenames)==2:
            flex_filename = outputfilenames[0]
            rigid_filename = outputfilenames[1]
        else:
            if self.debug: print "using default filenames: ", 
            flex_filename = self.ligand.name+'_'+self.molecule.name+'_flex.pdbqt'
            rigid_filename = self.ligand.name+'_'+self.molecule.name+'_rigid.pdbqt'
            if self.debug: print flex_filename, ' and ', rigid_filename
        self.write_flex(self.ligand, self.flex_residues, flex_filename)
        self.write_rigid(self.molecule, rigid_filename)
        delattr(self.flex_residues.atoms, 'used')
        delattr(self.flex_residues, 'setup')


    def write_rigid(self, molecule, outputfilename):
        #check each atom, write if not previously written
        outptr = open(outputfilename, 'w')
        for at in molecule.chains.residues.atoms:
            if not hasattr(at, 'used') or at.used==0:
                self.writer.write_atom(outptr, at)
        outptr.close()


    #support either ligand, a MolKit molecule or ligandfile, pdbqt file
    def write_flex(self, flex_residues, ligandfile, outputfilename):
        #need to know total rotatable bonds in flex_residues
        res_total = 0
        for r in flex_residues:
            res_total = res_total + r.torscount
        outfileptr = open(outputfilename,'w')
        ligfileptr = open(ligandfile)
        #this counter is used to reset the atom numbers
        # in the flexible residues and must be 1-based
        self.outatom_counter = 1
        for line in ligfileptr:
            #don't write any CONECT records
            if string.find(line, 'active torsions')>-1:
                lig_torscount = int(line.split()[1])
                self.torsionCtr = lig_torscount
                total = lig_torscount + res_total
                line = 'REMARK  %2d active torsions:\n' %(total)
                #print 'activeTorsionsLine= ', item
            if string.find(line, "ATOM")>-1:
                self.outatom_counter += 1
            if string.find(line, "HETA")>-1:
                self.outatom_counter += 1
            if string.find(line, 'CONECT')>-1: continue
            if string.find(line, 'MASTER')>-1: continue
            outfileptr.write(line)
        ligfileptr.close()
        for res in flex_residues:
            self.writeResidue(res, outfileptr)
        outfileptr.close()
        if self.debug: print "wrote flex file", outputfilename


    def writeResidue(self, res, outfileptr):
        # should have ALREADY eliminated residues with no active torsions
        # need to have already done autotors stuff to this residue's sidechain:
        # atoms should know which charge to use
        # also must have these fields: torscount, torsdof, and bonds know if active
        # also res.bondlist
        #start output
        outfileptr.write("BEGIN_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))
        #first write out remarks about torsions
        outfileptr.write("REMARK  " + "%d" %(res.torscount) + " active torsions:\n")
        outfileptr.write("REMARK  status: ('A' for Active; 'I' for Inactive)\n")
        #only want to process bonds pertaining to sideChain atoms
        for b in res.bondlist:
            if b.activeTors == 1: bTors = 'A'
            else:   bTors = 'I'
            if b.possibleTors:
                if bTors=='A':
                    self.torsionCtr = self.torsionCtr + 1
                    outstring = "REMARK " +" %3d  %s    between atoms: %-3s  and  %-3s \n" %(self.torsionCtr,bTors,b.atom1.name,b.atom2.name)
                else:
                    outstring = "REMARK " +"      %s    between atoms: %-3s  and  %-3s \n" %(bTors,b.atom1.name,b.atom2.name)
                outfileptr.write(outstring)
        #next write out  root which is always the CA atom
        outfileptr.write("ROOT\n")
        assert hasattr(res, 'rootlist')
        #reset used field to serve as visited flag
        res.atoms.used = 0
        #rootlist grows to include atoms up to first active tors in each subtree
        for at in res.rootlist:
            at.used = 1
            for bond in at.bonds:
                if bond.activeTors and bond.possibleTors: continue
                at2 = bond.atom1
                if at2==at: at2 = bond.atom2
                #only track and write out bonds to the sideChain atoms
                if at2 not in res.sideChain: continue
                if at2.used: continue
                if at2 not in res.rootlist:
                    res.rootlist.append(at2)
                    at2.rnum = len(res.rootlist)
        #remove all atoms which have no rotatable bonds
        #from BOTH the rootlist and the sideChain atoms
        #this means they will be written in the rigid portion
        badList = AtomSet([])
        for at in res.rootlist:
            hasTorsion=0
            for b in at.bonds:
                if b.activeTors:
                    hasTorsion=1
                    break
            if not hasTorsion:
                badList.append(at)
        if len(badList) and len(badList)!=len(res.rootlist):
            res.rootlist = res.rootlist - badList
            res.sideChain = res.sideChain - badList

        #now visit atoms connect to expanded rootlist
        for at in res.rootlist:
            at.number = self.outatom_counter
            self.writer.write_atom(outfileptr, at)
            at.newindex = self.outatom_counter
            at.used = 1
            self.outatom_counter = self.outatom_counter + 1
        outfileptr.write("ENDROOT\n")
        #last write out the rest of the stuff, using writeSubtree.....
        for at in res.rootlist:
            for b in at.bonds:
                at2 = b.atom1
                if at2 == at: at2 = b.atom2
                if at2 in res.rootlist:
                    continue
                if at2 not in res.sideChain:
                    at2.used = 0
                    continue
                if at2.used:
                    continue
                marker = self.outatom_counter
                outstring = "BRANCH %3d %3d\n"%(at.newindex ,marker)
                outfileptr.write(outstring)
                self.writeSubtree(outfileptr, at, at2)
                outstring = "ENDBRANCH %3d %3d\n"%(at.newindex,marker)
                outfileptr.write(outstring)
        outfileptr.write("END_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))



    def writeSubtree(self,outfptr, fromAtom, startAtom):
        """
            None <-writeSubtree(outfptr, fromAtom, startAtom)
            writeSubtree recursively visits the atoms in the current 
            'subtree' of the molecule in a Depth First Order traversal. 
            It is used to write out the molecule with the correct format 
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH 
            statements are added. 
        """
        if startAtom.used==0:
            startAtom.used = 1
            charge = startAtom.charge
            at = startAtom
            at.newindex = self.outatom_counter
            at.number=self.outatom_counter
            self.writer.write_atom(outfptr,at)
            self.outatom_counter = self.outatom_counter + 1
            marker = self.outatom_counter
            #outfptr.write(outstring)
            for bond in startAtom.bonds:
                nextAtom = bond.atom1
                if nextAtom == startAtom: nextAtom = bond.atom2
                if nextAtom == fromAtom: continue
                elif not nextAtom.used:
                    if bond.activeTors:
                        outstring = "BRANCH %3d %3d\n"%(marker-1,marker)
                        outfptr.write(outstring)
                    self.writeSubtree(outfptr, startAtom, nextAtom)
                    if bond.activeTors:
                        outstring = "ENDBRANCH %3d %3d\n"%(marker-1,marker)
                        outfptr.write(outstring)
        return

        
    def process_residues(self, residues):
        processed_residues = []
        for res in residues:
            if res.type=="PRO":
                continue
            if res.type=="HOH":
                continue
            ntors = self.setAutoFlexFields(res)
            #only keep residues with at least 1 torsion
            if ntors>0:
                processed_residues.append(res)
        return ResidueSet(processed_residues)
        

    def setAutoFlexFields(self, res):
        if self.debug: print "in setAutoFlexFields with ", res.full_name()
        if hasattr(res, 'setup'): 
            return
        res.setup = 1
        res.atoms.used = 0
        res.atoms.bonds[0].possibleTors = 0
        res.atoms.bonds[0].activeTors = 0
        #backbone_names = ['C','N','O','HN','HN1','HN2', 'HA', 
        #            'H1','H2','H3','HO', 'H']
        #sidechain = res.atoms - res.atoms.get('backbone+h')
        #restore CA
        #sidechain = res.sideChain = sidechain + res.atoms.get('CA')
        #sidechain = res.atoms.get(lambda x: x.name not in backbone_names)
        sidechain = res.atoms.get('sidechain')
        ca_atom = res.atoms.get('CA')
        if ca_atom is not None:
            sidechain = sidechain + ca_atom
            
        res.sideChain = sidechain
        bondlist = res.bondlist = sidechain.bonds[0]
        bondlist.leaf = 0
        bondlist.possibleTors = 0
        bondlist.activeTors = 0
        rbs = RotatableBondSelector()
        rotatables = rbs.select(bondlist)
        if self.debug: print "len(rotatables)=", len(rotatables)
        for b in rotatables:
            b.possibleTors = 1
            b.activeTors = 1
        #turn off amide and guanidinium unless specifically allowed
        if 'amide' not in self.allowed_bonds:
            amide_bnds = AmideBondSelector().select(rotatables)
            for b in amide_bnds:
                b.possibleTors = 1
                b.activeTors = 0
        if 'guanidinium' not in self.allowed_bonds:
            guanidinium_bnds = GuanidiniumBondSelector().select(rotatables)
            for b in guanidinium_bnds:
                b.possibleTors = 1
                b.activeTors = 0
        for b in self.non_rotatable_bonds:
            b.activeTors = 0
        res.torscount = len(rotatables.get(lambda x:x.activeTors==1))
        #this field is not used in AutoDock4
        res.torsdof = res.torscount
        if ca_atom is not None:
            res.rootlist = ca_atom
            res.root = ca_atom[0]
        else:
            at0 = res.atoms.get(lambda x: x._uniqIndex==0)[0]
            res.rootlist = AtomSet([at0])
            res.root = at0
        return res.torscount



class AD4FlexibleReceptorPreparation:
    """
    Facade for preparing 'receptor' input for an AutoDock4 experiment 
    with flexible sidechains in the receptor in flexres_filename
    rest of the receptor in the rigid_filename

    This preparation involves:
        1. adding autodock_element types 
        2. selecting flexible residues 
        3. specifying patterns of flexiblilty in side-chains of selected
        residues
        4. writing flexres file containing formatted residues
        5. writing a rigid pdbqt outputfile containing remaining receptor
        atoms


    Receptor preparation collaborators extend baseclass collaborators:
        in AutoDockTools:
            in atomTypeManager classes for conforming to AD atomtypes:
                AutoDock4_AtomTyper
       in this file:
            AD4ReceptorWriter
        
"""

    def __init__(self, receptor, mode='automatic',
                    rigid_filename=None,  residues=[],
                    allowed_bonds='',  #backbone on, amide + guanidinium off
                    non_rotatable_bonds=[], flexres_filename=None, 
                    reassign=True, debug=False):

        #?NEED TO TYPE ATOMS?
        self.receptor = receptor
        self.debug = debug
        file_type = os.path.splitext(os.path.basename(receptor.parser.filename))[1]
        if file_type=='.pdbqt':
            pass
        elif file_type in ['.pdbqs','.pdbq','.mol2', '.pdb']:
            if file_type in ['.pdbqs','.pdbq']:
                try:
                    delattr(receptor.allAtoms, 'autodock_element')
                except:
                    pass
            atomTyper = AutoDock4_AtomTyper()
            atomTyper.setAutoDockElements(receptor, reassign=reassign) #
        else:
            print file_type, " unknown type: must be .pdbqt, .pdbqs, .pdbq, .pdb or .mol2" 
        self.torsion_ct = 0
        self.writer = PdbqtWriter()
        #process the flexres
        try:
            self.allowed_bonds = allowed_bonds.split('_')  # a string "backbone"
        except:
            self.allowed_bonds = ['']
        self.non_rotatable_bonds = non_rotatable_bonds
        if self.debug: print "about to process_residues: ", residues
        flex_residues = self.flex_residues = self.process_residues(residues)
        #MODE SWITCH 5: write outputfiles     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: are there other records which should be written?
            if flexres_filename is None:
                flexres_filename = self.receptor.name+'_flex.pdbqt'
            if rigid_filename is None:
                rigid_filename = self.receptor.name+'_rigid.pdbqt'
            self.write_flex(flex_residues, flexres_filename)
            if self.debug: print "wrote flexible file:", flexres_filename
            self.write_rigid(self.receptor, rigid_filename)
            if self.debug: print "wrote rigid file:", rigid_filename


    def write(self, outputfilenames=None):
        if outputfilenames is not None and len(outputfilenames)==2:
            flexres_filename = outputfilenames[0]
            rigid_filename = outputfilenames[1]
        else:
            if self.debug: print "using default filenames: ", 
            flexres_filename = self.receptor.name+'_flex.pdbqt'
            rigid_filename = self.receptor.name+'_rigid.pdbqt'
            if self.debug: print flexres_filename, ' and ', rigid_filename
        self.write_flex(self.flex_residues, flexres_filename)
        self.write_rigid(self.receptor, rigid_filename)
        delattr(self.flex_residues.atoms, 'used')
        delattr(self.flex_residues, 'setup') 

    def write_rigid(self, molecule, rigid_filename):
        #check each atom, write if not previously written
        outptr = open(rigid_filename, 'w')
        for at in molecule.chains.residues.atoms:
            if not hasattr(at, 'used') or at.used==0:
                self.writer.write_atom(outptr, at)
        outptr.close()


    def write_flex(self, flex_residues, flexres_filename):
        #need to know total rotatable bonds in flex_residues
        if self.debug: 
            print "in write_flex: len(flex_residues)=", len(flex_residues), 'ffn=', flexres_filename
        res_total = 0
        for r in flex_residues:
            res_total = res_total + r.torscount
        if self.debug: print "res_total=", res_total
        self.writtenAtoms = []
        outfileptr = open(flexres_filename,'w')
        self.outatom_counter = 1
        self.outatom_counter = 1
        self.torsionCtr = 0
        for item in flex_residues:    
            self.writeResidue(item, outfileptr)
        outfileptr.close()
        if self.debug: print "wrote flex file", flexres_filename


    def writeResidue(self, res, outfileptr):
        # should have ALREADY eliminated residues with no active torsions
        # need to have already done autotors stuff to this residue's sidechain:
        # atoms should know which charge to use
        # also must have these fields: torscount, torsdof, and bonds know if active
        # also res.bondlist
        #start output
        if self.debug: print "in writeResidue with ", res.full_name()
        outfileptr.write("BEGIN_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))
        #first write out remarks about torsions
        outfileptr.write("REMARK  " + "%d" %(res.torscount) + " active torsions:\n")
        outfileptr.write("REMARK  status: ('A' for Active; 'I' for Inactive)\n")
        #only want to process bonds pertaining to sideChain atoms
        for b in res.bondlist:
            if b.activeTors == 1: bTors = 'A'
            else:   bTors = 'I'
            if b.possibleTors:
                if bTors=='A':
                    self.torsion_ct += 1
                    outstring = "REMARK " +" %3d  %s    between atoms: %-3s  and  %-3s \n" %(self.torsion_ct,bTors,b.atom1.name,b.atom2.name)
                else:
                    outstring = "REMARK " +"      %s    between atoms: %-3s  and  %-3s \n" %(bTors,b.atom1.name,b.atom2.name)
                outfileptr.write(outstring)
        #next write out  root which is always the CA atom
        outfileptr.write("ROOT\n")
        assert hasattr(res, 'rootlist')
        #reset used field to serve as visited flag
        res.atoms.used = 0
        #rootlist grows to include atoms up to first active tors in each subtree
        for at in res.rootlist:
            at.used = 1
            for bond in at.bonds:
                if bond.activeTors and bond.possibleTors: continue
                at2 = bond.atom1
                if at2==at: at2 = bond.atom2
                #only track and write out bonds to the sideChain atoms
                if at2 not in res.sideChain: continue
                if at2.used: continue
                if at2 not in res.rootlist:
                    res.rootlist.append(at2)
                    at2.rnum = len(res.rootlist)
        #remove all atoms which have no rotatable bonds
        #from BOTH the rootlist and the sideChain atoms
        #this means they will be written in the rigid portion
        badList = AtomSet([])
        for at in res.rootlist:
            hasTorsion=0
            for b in at.bonds:
                if b.activeTors:
                    hasTorsion=1
                    break
            if not hasTorsion:
                badList.append(at)
        if len(badList) and len(badList)!=len(res.rootlist):
            res.rootlist = res.rootlist - badList
            res.sideChain = res.sideChain - badList
        #now visit atoms connect to expanded rootlist
        for at in res.rootlist:
            at.number = self.outatom_counter
            self.writer.write_atom(outfileptr, at)
            at.newindex = self.outatom_counter
            at.used = 1
            self.outatom_counter = self.outatom_counter + 1
        outfileptr.write("ENDROOT\n")
        #last write out the rest of the stuff, using writeSubtree.....
        for at in res.rootlist:
            for b in at.bonds:
                at2 = b.atom1
                if at2 == at: at2 = b.atom2
                if at2 in res.rootlist:
                    continue
                if at2 not in res.sideChain:
                    continue
                if at2.used:
                    continue
                self.process(at, at2, outfileptr)
                #marker = self.outatom_counter
                #outstring = "BRANCH %3d %3d\n"%(at.newindex ,marker)
                #outfileptr.write(outstring)
                #self.writeSubtree(outfileptr, at, at2)
                #outstring = "ENDBRANCH %3d %3d\n"%(at.newindex,marker)
                #outfileptr.write(outstring)
        outfileptr.write("END_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))


    def process(self, fromAtom, nextAtom, outfileptr):
        startIndex = fromAtom.number
        endIndex = self.outatom_counter
        outstring = "BRANCH %3d %3d\n"%(startIndex, endIndex)
        outfileptr.write(outstring)
        queue = self.writeBreadthFirst(outfileptr, fromAtom, nextAtom)
        if self.debug: print fromAtom.name, ':', nextAtom.name, ': queue=', queue
        if len(queue):
            for fromAtom, nextAtom in queue:
                if self.debug: print " processing queue entry: ", fromAtom.name, '-', nextAtom.name
                self.process(fromAtom, nextAtom, outfileptr)
        outstring = "ENDBRANCH %3d %3d\n"%(startIndex, endIndex)
        outfileptr.write(outstring)

    def writeLevel(self, atom, outfptr):
        """
        write all atoms bonded to atoms bonded to this atom by non-rotatable
        bonds
        """
        if self.debug: 
            print "\n\nin writeLevel with ", atom.name, " outatom_counter=", self.outatom_counter
            print "len(", atom.name, ").bonds=", len(atom.bonds)
        queue = []
        nextAts = []
        for b in atom.bonds:
            if self.debug:
                print "processing b=", b.atom1.name, '-', b.atom2.name, ' activeTors=', b.activeTors
                print 'atom1 in writtenAtoms=', b.atom1 in self.writtenAtoms
                print 'atom2 in writtenAtoms=', b.atom2 in self.writtenAtoms
            if b.activeTors: 
                at2 = b.atom1
                if at2==atom: at2=b.atom2
                queue.append((atom, at2))
                if self.debug: print atom.name, 'wL: queue=', queue
                continue
            a2 = b.atom1
            if a2==atom:
                a2 = b.atom2
            if a2.used: 
                if self.debug: print "!!a2 is already used!!", a2.name
                continue    
            if a2 not in self.writtenAtoms:
                a2.number = a2.newindex = self.outatom_counter
                if self.debug: print "writeLevel: wrote bonded atom named=", a2.name, 'a2.used=', a2.used
                self.writer.write_atom(outfptr, a2)
                self.writtenAtoms.append(a2)
                a2.used = 1
                self.outatom_counter+=1
                nextAts.append(a2)
        for a2 in nextAts:
            if self.debug: 
                print 'in for nextAts loop with a2=', a2.name
                print 'calling wL'
            nq = self.writeLevel(a2, outfptr)
            if len(nq):
                if self.debug: print "extending queue with", nq
                queue.extend(nq)
        if self.debug:
            print 'returning queue=', queue
        return queue
            

    def writeBreadthFirst(self, outfptr, fromAtom, startAtom):
        """
            None <-writeBreadthFirst(outfptr, fromAtom, startAtom)
            writeBreadthFirst visits all the atoms in the current level
            then the first level down etc in a Breadth First Order traversal. 
                            1                <-1
                        5 6   7 8            <-3
                     9 10   11 12            <-4
            It is used to write out the molecule with the correct format 
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH 
            statements are added. 
        """
        if self.debug: 
            print "in wBF with fromAtom=", fromAtom.name, '+ startAtom=', startAtom.name, 'startAtom.used=', startAtom.used
        queue = []
        if startAtom.used==0:
            startAtom.used = 1
            startAtom.number = startAtom.newindex = self.outatom_counter
            self.writer.write_atom(outfptr,startAtom)
            if self.debug: print 'wBF: wrote ', startAtom.name
            self.writtenAtoms.append(startAtom)
            self.outatom_counter += 1
            if self.debug: print "self.outatom_counter=", self.outatom_counter
            activeTors = []
            #outfptr.write(outstring)
            for bond in startAtom.bonds:
                if not hasattr(bond, 'activeTors'):
                    continue
                at2 = bond.atom1
                if at2==startAtom: at2 = bond.atom2
                if at2==fromAtom: continue  #this is current bond
                elif not at2.used:
                    if bond.activeTors:
                        queue.append((startAtom,at2))
                    else:
                        at2.number = at2.newindex = self.outatom_counter
                        if self.debug: 
                            print "\n\nwriting and calling wL with nA=", at2.name, '-', at2.number
                        self.writer.write_atom(outfptr, at2)
                        if self.debug: print 'wBF2: wrote ', at2.name
                        at2.written = 1
                        self.writtenAtoms.append(at2)
                        at2.newindex = self.outatom_counter
                        self.outatom_counter = self.outatom_counter + 1
                        if self.debug: print '!!!2:calling wL'
                        newQ = self.writeLevel(at2, outfptr)
                        if self.debug: print "newQ=", newQ
                        at2.used = 1
                        if len(newQ): 
                            if self.debug: print "@@@@len(newq)=", len(newQ)
                            queue.extend(newQ)
                            if self.debug: print "queue=", queue
            if self.debug: 
                print " currently queue=",
                for atom1, atom2 in queue: 
                    print atom1.name, '-', atom2.name, ',',
                print
        return  queue


    def writeSubtree(self,outfptr, fromAtom, startAtom):
        """
            None <-writeSubtree(outfptr, fromAtom, startAtom)
            writeSubtree recursively visits the atoms in the current 
            'subtree' of the molecule in a Depth First Order traversal. 
            It is used to write out the molecule with the correct format 
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH 
            statements are added. 
        """
        if startAtom.used==0:
            startAtom.used = 1
            at = startAtom
            for bond in startAtom.bonds:
                if bond.activeTors: 
                    continue
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom: 
                    nextAtom = bond.atom2
                if nextAtom==fromAtom: 
                    continue
                if not nextAtom.used:
                    if hasattr(bond,'incycle'):
                        if not hasattr(nextAtom, 'cycleout'):
                            nextAtom.cycleout = 1
                            nextAtom.newindex = self.outatom_counter
                            nextAtom.number = self.outatom_counter
                            self.writer.write_atom(outfptr,nextAtom)
                            self.writtenAtoms.append(nextAtom)
                            self.outatom_counter = self.outatom_counter+1
                    else:
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
            for bond in startAtom.bonds:
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom: 
                    nextAtom = bond.atom2
                if nextAtom==fromAtom: 
                    continue
                if not nextAtom.used:
                    testcond = len(nextAtom.bonds)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "BRANCH %3d %3d\n"%(at.newindex,marker)
                            outfptr.write(outstring)
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
                    self.WriteSubtree(startAtom, nextAtom)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "ENDBRANCH %3d %3d\n"%(at.newindex,marker)
                            outfptr.write(outstring)
        return


#        if startAtom.used==0:
#            startAtom.used = 1
#            charge = startAtom.charge
#            at = startAtom
#            at.newindex = self.outatom_counter
#            at.number=self.outatom_counter
#            self.PdbqtWriter.write_atom(outfptr,at)
#            self.outatom_counter = self.outatom_counter + 1
#            marker = self.outatom_counter
#            #outfptr.write(outstring)
#            for bond in startAtom.bonds:
#                nextAtom = bond.atom1
#                if nextAtom == startAtom: nextAtom = bond.atom2
#                if nextAtom == fromAtom: continue
#                elif not nextAtom.used:
#                    if bond.activeTors:
#                        outstring = "BRANCH %3d %3d\n"%(marker-1,marker)
#                        outfptr.write(outstring)
#                    self.writeSubtree(outfptr, startAtom, nextAtom)
#                    if bond.activeTors:
#                        outstring = "ENDBRANCH %3d %3d\n"%(marker-1,marker)
#                        outfptr.write(outstring)
#        return

        
    def process_residues(self, residues):
        processed_residues = []
        for res in residues:
            if res.type=="PRO":
                continue
            if res.type=="HOH":
                continue
            if self.debug: print 'calling setAutoFlexFields with ', res.full_name()
            ntors = self.setAutoFlexFields(res)
            if self.debug:  print "ntors=", ntors
            #only keep residues with at least 1 torsion
            if ntors>0:
                processed_residues.append(res)
        return ResidueSet(processed_residues)
        

    def setAutoFlexFields(self, res):
        if self.debug: print "in setAutoFlexFields with ", res.full_name()
        if hasattr(res, 'setup'): 
            return
        res.setup = 1
        res.atoms.used = 0
        res.atoms.bonds[0].possibleTors = 0
        res.atoms.bonds[0].activeTors = 0
        #backbone_names = ['C','N','O','HN','HN1','HN2', 'HA', 
        #            'H1','H2','H3','HO', 'H']
        #sidechain = res.atoms - res.atoms.get('backbone+h')
        #restore CA
        #sidechain = res.sideChain = sidechain + res.atoms.get('CA')
        #sidechain = res.atoms.get(lambda x: x.name not in backbone_names)
        sidechain = res.atoms.get('sidechain')
        ca_atom = res.atoms.get('CA')
        if ca_atom is not None:
            sidechain = sidechain + ca_atom
        res.sideChain = sidechain
        bondlist = res.bondlist = sidechain.bonds[0]
        bondlist.leaf = 0
        bondlist.possibleTors = 0
        bondlist.activeTors = 0
        rbs = RotatableBondSelector()
        rotatables = rbs.select(bondlist)
        if self.debug: print "len(rotatables)=", len(rotatables)
        for b in rotatables:
            b.possibleTors = 1
            b.activeTors = 1
        #turn off amide and guanidinium unless specifically allowed
        if 'amide' not in self.allowed_bonds:
            amide_bnds = AmideBondSelector().select(rotatables)
            for b in amide_bnds:
                b.possibleTors = 1
                b.activeTors = 0
        if 'guanidinium' not in self.allowed_bonds:
            guanidinium_bnds = GuanidiniumBondSelector().select(rotatables)
            for b in guanidinium_bnds:
                b.possibleTors = 1
                b.activeTors = 0
        for b in self.non_rotatable_bonds:
            b.activeTors = 0
        res.torscount = len(rotatables.get(lambda x:x.activeTors==1))
        #this field is not used in AutoDock4
        res.torsdof = res.torscount
        if ca_atom is not None:
            res.rootlist = ca_atom
            res.root = ca_atom[0]
        else:
            at0 = res.atoms.get(lambda x: x._uniqIndex==0)[0]
            res.rootlist = AtomSet([at0])
            res.root = at0
        if self.debug: 
            print 'returning res.torscount=', res.torscount
            print 'res.root=', res.root.name
        return res.torscount

                
if __name__ == '__main__':
    from MolKit import Read
    m = Read("/mgl/work4/rhuey/dev23/1hsg.pdb")[0]
    #m = Read("/mgl/work4/rhuey/dev23/MolKit/Tests/Data/protease.pdb")[0]
    m.buildBondsByDistance()
    import os
    #assert os.path.exists("protease.pdbqs")==False
    RPO = ReceptorPreparation(m, debug=True)
    #assert os.path.exists("protease.pdbqs")==False
    os.system('date')
    RPO.write()
    os.system('date')
    #assert os.path.exists("protease.pdbqs")==True


#    m2 = Read("/mgl/work4/rhuey/dev23/MolKit/Tests/Data/indinavir.pdb")[0]
#    m2.buildBondsByDistance()
#    LPO = LigandPreparation(m2)
#    import os
#    #assert os.path.exists("indinavir.pdbq")==False
#    LPO.write()
#    #assert os.path.exists("indinavir.pdbq")==True
