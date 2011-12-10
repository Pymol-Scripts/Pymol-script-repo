from src.forcefield import *
from peoe_PDB2PQR import PEOE as calc_charges
from src.pdb import *
from src.definitions import *
from pdb2pka import NEWligand_topology
import string

def initialize(definition, ligdesc, pdblist, verbose=0):
    """
        Initialize a ligand calculation by adding ligand atoms to the definition.
        This code adapted from pdb2pka/pka.py to work with existing PDB2PQR code.

        Parameters
            definition:  The definition object for PDB2PQR
            ligdesc:    A file desciptor the desired ligand file (string)
            pdblist:    A list of objects from the PDB module in the original protein (list)
            verbose:    When 1, script will print information to stdout

        Returns
            protein:   The protein created by the original atoms and the ligand
            definition: The updated definition for this protein
            Lig:       The ligand_charge_handler object
    """
    
    Lig = ligand_charge_handler()
    Lig.read(ligdesc)
    
    # Create the ligand definition from the mol2 data

    MOL2FLAG = True
    X=NEWligand_topology.get_ligand_topology(Lig.lAtoms,True)

    # Add it to the 'official' definition

    ligresidue = DefinitionResidue()
    ligresidue.name = "LIG"
    i = 1
    atommap = {}
    for line in X.lines[:-2]:
        obj = DefinitionAtom()
        entries = string.split(line)
        if len(entries) != 4:
            raise ValueError, "Invalid line for MOL2 definition!"
        name = entries[0]
        obj.name = name
        obj.x = float(entries[1])
        obj.y = float(entries[2])
        obj.z = float(entries[3])
        atommap[i] = name
        ligresidue.map[name] = obj
        i += 1

    # The second to last line has a list of bonded partners

    line = X.lines[-2]
    bonds = string.split(line)

    for i in range(0,len(bonds),2):
        bondA = int(bonds[i])
        bondB = int(bonds[i+1])
        atomA = ligresidue.getAtom(atommap[bondA])
        atomB = ligresidue.getAtom(atommap[bondB])

        atomA.bonds.append(atommap[bondB])
        atomB.bonds.append(atommap[bondA])

    # The last line is not yet supported - dihedrals
    
    definition.map["LIG"] = ligresidue

    # Look for titratable groups in the ligand

    ligand_titratable_groups=X.find_titratable_groups()

    if verbose:
        print "ligand_titratable_groups", ligand_titratable_groups
    
    # Append the ligand data to the end of the PDB data

    newpdblist=[]
       
    # First the protein
    nummodels = 0 
    for line in pdblist:
        if isinstance(line, END) or isinstance(line,MASTER): continue
        elif isinstance(line, MODEL):
            nummodels += 1
            if nummodels > 1: break
        else:
            newpdblist.append(line)
    
    # Now the ligand

    for e in Lig.lAtoms: 
        e.resName = "LIG"
        newpdblist.append(e)
     
    protein = Protein(newpdblist, definition)
    for rrres in  protein.chainmap['L'].residues:
        for aaat in rrres.atoms:
            for ligatoms in Lig.lAtoms:
                if ligatoms.name == aaat.name:
                    aaat.sybylType = ligatoms.sybylType
    
                    # setting the formal charges

                    if ligatoms.sybylType == "O.co2":
                        aaat.formalcharge = -0.5
                    else: aaat.formalcharge = 0.0
                    xxxlll = []
                    #for xxx in ligatoms.lBondedAtoms:
                    for bond in ligresidue.getAtom(aaat.name).bonds:
                        xxxlll.append(bond)
        
                    aaat.intrabonds = xxxlll
   
                    # charge initialisation must happen somewhere else
                    aaat.charge = 0.0

                    #print aaat.name, aaat.charge, aaat.formalcharge, aaat.intrabonds

    return protein, definition, Lig


class ligforcefield(Forcefield):
    """
    Derived class of  forcefield.py for charge assignment on ligands
    """

    def __init__(self, ff, lig_instance):
        """
            Initialize the class by parsing the definition file

            Parameters
                ff: The name of the forcefield (string)

            Additionally, ligands can be considered within this class
        """
        print "lig_instance", lig_instance
        self.residues = {}
        self.name = ff
        defpath = ""
        if ff == "amber":
            defpath = AMBER_FILE
        elif ff == "charmm":
            defpath = CHARMM_FILE
        elif ff == "parse":
            defpath = PARSE_FILE
        else:
            raise ValueError, "Invalid forcefield %s!" % ff

        if not os.path.isfile(defpath):
            for path in sys.path:
                testpath = "%s/%s" % (path, defpath)
                if os.path.isfile(testpath):
                    defpath = testpath
                    break
        if not os.path.isfile(defpath):
            raise ValueError, "Unable to find forcefield %s!" % defpath

        file = open(defpath)
        lines = file.readlines()
        for line in lines:
            if not line.startswith("#"):
                fields = string.split(line)
                resname = fields[0]
                atomname = fields[1]
                charge = float(fields[2])
                radius = float(fields[3])

                atom = ForcefieldAtom(atomname, charge, radius)

                myResidue = self.getResidue(resname)
                if myResidue == None:
                    myResidue = ForcefieldResidue(resname)
                    self.residues[resname] = myResidue
                myResidue.addAtom(atom)
            #
        ### PC - charge assignment on ligand
        ###
        #self.lig = MOL2MOLECULE()
        self.lig=lig_instance #lig_shit()
        #self.lig.read(ligfilename)
        return

    
    def getParams(self,residue,name):
        """
            Get the parameters associated with the input fields.
            The residue itself is needed instead of simply its name
            because  the forcefield may use a different residue name
            than the standard amino acid name.

            Parameters
                residue:  The residue (residue)
                name:     The atom name (string)
            Returns
                charge:   The charge on the atom (float)
                radius:   The radius of the atom (float)
        """
        charge = None
        radius = None
        resname = ""
        atomname = ""
        #
        # PC - 230306 - we need to put the setting of formal charges in another place
        #for at in self.lig.lAtoms:
        #    at.charge = 0.0
        if self.name == "amber" and residue.type != 2:
            resname, atomname = self.getAmberParams(residue, name)
        elif self.name == "charmm" and residue.type != 2:
            resname, atomname = self.getCharmmParams(residue, name)
        elif self.name == "parse" and residue.type != 2:
            resname, atomname = self.getParseParams(residue, name)
        defresidue = self.getResidue(resname)
        ### This is a rather quick and dirty solution
        if residue.type == 2:
            charge,radius = self.getChargeAndRadius(residue,name)
        if defresidue != None:
            atom = defresidue.getAtom(atomname)
        else:
            atom = None
        if atom != None:
            #print "XXX___LigandAtom___XXX", atom.name # PC 050506
            charge = atom.get("charge")
            radius = atom.get("radius")
        return charge, radius

    def getChargeAndRadius(self,residue,atomname):
        self.lig.make_up2date(residue)
        return self.lig.ligand_props[atomname]['charge'],self.lig.ligand_props[atomname]['radius']


#
# ---------
#

BondiiRadiiDict = {"C": 1.70,
                   "N": 1.50,
                   "O": 1.40,
                   "S": 1.85,
                   "H": 1.05,
                   "Br":2.50,
                   "F": 1.20,
                   "P": 1.90}

class ligand_charge_handler(MOL2MOLECULE):
    """Make sure that we are up to date with respect to the charge calculation"""

    def make_up2date(self,residue):
        #
        # Check if the structure of the ligand is
        # identical to the one we have
        #
        if not getattr(self,'ligand_props',None):
           self.recalc_charges(residue)
           qqqgesges = 0.0
           for aa in residue.atoms:
               #print "newly_calced  %s  %1.4f " %(aa.name, aa.charge)
               qqqgesges = qqqgesges +  aa.charge
           #print "-------------------------------"
           #print "newly_calced - net charge %1.4f" %(qqqgesges)
           #print
        else:
            atoms_last_calc=self.ligand_props.keys()
            #
            # Get the atoms presently in the pdb2pqr array
            #
            atoms_now=[]
            for at in residue.atoms:
                atoms_now.append(at.name)
            atoms_now.sort()
            #print "atoms_now ",atoms_now
            atoms_last_calc.sort()
            #
            xxnetqxx = 0.0
            for aa in residue.atoms:
                #print "NOT recalced  %s  %1.4f  %1.4f " %(aa.name, aa.charge, aa.formalcharge)
                xxnetqxx = xxnetqxx+aa.charge
            #print "NOT recalced, net_q : %1.2f" %(xxnetqxx)
            #print "###########################"
            #
            # Did we calculate charges for exactly this set of atoms?
            #
            #print "atoms_last_calc        ", atoms_last_calc
            #print "Atoms_this_time_around ", atoms_now
            if atoms_now!=atoms_last_calc:
                #
                # No - Recalc charges
                #
                for atom in atoms_now:
                    if not atom in atoms_last_calc:
                        print 'This atom was missing before',atom
                        print 'If it is a hydrogen for the titratable, we need to create a bond entry!'
                        # We should be here only if is a titratable
                        for current_atom in atoms_now:
                            # check if it't a titratable H
                            for res_atoms in residue.atoms:
                                if current_atom == res_atoms.name and "titratableH"  in dir(res_atoms):
                                    print "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
                                    print "been here"
                                    for ResAtoms in residue.atoms:
                                        ResAtoms.formalcharge = 0.0
                                    self.recalc_charges(residue)
                for atom in atoms_last_calc:
                    if not atom in atoms_now:
                        print 'This atom used to be here, but is now missing',atom
                #self.recalc_charges(residue)
                xxxnetq = 0.0
                for xxx in residue.atoms:
                    print "after neutralizing %s  %1.4f" %(xxx.name, xxx.charge)
                    xxxnetq = xxxnetq+xxx.charge
                print '-----------------------'
                print "net charge: %1.4f" % (xxxnetq)
                print
                print
            else:
                # Yes - nothing to do
                pass
        #
        # Now we have the charges
        #
        return

    #
    # ----
    #   

    def recalc_charges(self,residue):
        #
        # Recalculate the charges
        #
        self.ligand_props={}
        #
        # Here we have to update the atoms that are used for the charge calculation
        #
        # DO IT!!
        #
        #print "I am calling calc_charges"
        #
        # initialising the charges
        #
        for att in residue.atoms:
            att.charge= att.formalcharge
        #
        calc_charges(residue)
        #
        # Loop over all atoms
        #
        #print 'Atoms passed to Pauls routines'
        for at in residue.atoms: # WAS: self.lAtoms:
            ele = at.sybylType.split('.')[0]
            charge = at.charge
            if BondiiRadiiDict.has_key(ele):
                radius = BondiiRadiiDict[ele]
            else:
                raise 'radius not known for',ele
            #
            # Store the radii and charges
            #
            self.ligand_props[at.name]={'charge':charge,'radius':radius}
        return    
