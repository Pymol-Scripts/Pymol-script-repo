#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_covalent_flexres.py,v 1.7 2009/03/30 21:42:51 rhuey Exp $
#
import os 
from PyBabel.atomTypes import AtomHybridization
from mglutil.math.rigidFit import RigidfitBodyAligner
from MolKit import Read
from MolKit.molecule import AtomSet
from MolKit.bondSelector import RotatableBondSelector, AmideBondSelector
from MolKit.bondSelector import GuanidiniumBondSelector, LeafBondSelector
from mglutil.math.rotax import rotax
from math import sqrt

from AutoDockTools.MoleculePreparation import AD4FlexibleReceptorPreparation

try:
    import Numeric
except:
    import numpy.oldnumeric as Numeric

def dist(a1, a2):
    c1 = Numeric.array(a1.coords)
    c2 = Numeric.array(a2.coords)
    d = c2-c1
    return sqrt(sum(d*d))


def getBondedAtoms(atom):
    l = AtomSet()
    for b in atom.bonds:
        l.append(b.neighborAtom(atom))
    return l


def detect_clashes(atomset1, atomset2, cutoff=1.2, verbose=False):  #check this v. autodock/mdist.h
    #exclude bonded pairs
    clashes = []
    for at1 in atomset1:
        if verbose: print "processing ", at1.name,
        bonded = getBondedAtoms(at1)
        if verbose: print "len(bonded)=", len(bonded),
        for at2 in atomset2:
            if verbose: print "vs ", at2.name
            if at2 in bonded:
                continue
            d = dist(at1,at2)
            BondOrderSum = at1.bondOrderRadius + at2.bondOrderRadius
            if verbose: print 'distance = % 6.4f vs BondOrderSum*cutoff = % 6.4f' %(d, BondOrderSum*cutoff)
            if d <= BondOrderSum * cutoff:
                if verbose: print "%s - %s : % 6.4f - % 6.4f" %(at1.name, at2.name, d, BondOrderSum)
                clashes.append(((at1, at2),d))
    return clashes


def setAutoFlexFields(res):
    #process residues
    if hasattr(res, 'setup'): 
        return
    res.setup = 1
    res.atoms.used = 0
    res.atoms.bonds[0].possibleTors = 0
    res.atoms.bonds[0].activeTors = 0
    backbone_names = ['C','N','O','HN','HN1','HN2', 'HA', 
                'H1','H2','H3','HO', 'H']
    #includes CA
    sidechain = res.atoms.get(lambda x: x.name not in backbone_names)
    res.sideChain = sidechain
    bondlist = res.bondlist = sidechain.bonds[0]
    bondlist.leaf = 0
    bondlist.possibleTors = 0
    bondlist.activeTors = 0
    rbs = RotatableBondSelector()
    rotatables = rbs.select(bondlist)
    for b in rotatables:
        b.possibleTors = 1
        b.activeTors = 1
    amides = AmideBondSelector().select(bondlist)
    for b in amides:
        b.activeTors = 0
        b.possibleTors = 1
    guanidiniums = GuanidiniumBondSelector().select(bondlist)
    for b in guanidiniums:
        b.activeTors = 0
        b.possibleTors = 1
    leaves = LeafBondSelector().select(bondlist)
    for b in leaves:
        b.activeTors = 0
        b.possibleTors = 0
    res.torscount = len(bondlist.get(lambda x: x.activeTors==1))
    #this field is not used in AutoDock4
    res.torsdof = res.torscount
    res.torscount = len(bondlist.get(lambda x: x.activeTors==1))
    res.torsdof = res.torscount
    
    caAtoms = res.atoms.get(lambda x: x.name=='CA')
    #get returns an AtomSet
    if caAtoms:  #this checks for len(caAtoms)
        res.rootlist = caAtoms
    else:
        res.rootlist = AtomSet([res.atoms.get(lambda x: x._uniqIndex == 0)[0]])
        res.sideChain = res.atoms


def setRelativeTorsion(atom1, atom2, angle): #mov_atoms=None):
    #rat1, rat2, ??30degrees??
    mol = atom1.top
    if mov_atoms is None:
        mov_atoms = mol.subTree(atom1, atom2, mol.allAtoms)
    assert len(mov_atoms)
    mov_coords = Numeric.array(mov_atoms.coords)
    lenCoords = len(mov_coords)
    x = Numeric.array(atom1.coords)
    y = Numeric.array(atom2.coords)
    rot = (angle * 3.14159/180.)%(2 * Numeric.pi)
    matrix = rotax(x, y, rot)
    _ones = Numeric.ones(lenCoords, 'f')
    _ones.shape = (lenCoords,1)
    mov_coords = Numeric.concatenate((mov_coords, _ones),1)
    newcoords = Numeric.dot(mov_coords, matrix)
    nc = newcoords[:,:3].astype('f')
    for i in range(lenCoords):
        at = mov_atoms[i]
        at._coords[at.conformation] = nc[i].tolist()



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_covalent_flexres.py -r receptor_filename -l ligand_to_superimpose_filename -b atom1_atom2_atom3 -R residue_name"
        print "    Description of command..."
        print "         -r   receptor_filename"
        print "         -l   ligand_filename"
        print "         -b   ligand atom1_atom2_atom3 to superimpose (eg SD_CG_CB)"
        print "         -R   name of receptor residue onto which the bond atoms "
        print "                   in the ligand will be superimposed"
        print "    Optional parameters:"
        print "        [-B]  alternate atom1_atom2_atom3 in receptor residue to superimpose "
        print "            (default uses atom1_atom2_atom3 names specified with -b flag)"
        print "        [-v]  verbose output"
        print "        [-g pdbqt_filename] (rigid output filename)"
        print "        [-x pdbqt_filename] (flexible output filename)"
        print "                                              "
        print "  NOTE the order of the atoms specified for the bond is crucial:"   
        print "        '-b atom1_atom2_atom3' "
        print "   specifies this superposition:"
        print "       LIGAND-atom1-atom2-atom3"   
        print "              atom1-atom2-atom3-RECEPTOR"   
        print "  Process removes atom1-atom2-atom3-RECEPTOR bond and adds bond LIGAND-atom1-atom2-atom3-RECEPTOR"
        print "                                   ^                                                    ^"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:b:R:B:g:x:vh')
    except getopt.GetoptError, msg:
        print 'prepare_covalent_flexres.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptor
    receptor_filename =  None
    #-l: ligand
    ligand_filename =  None
    #-b: atom1_atom2_atom3
    atom1_atom2_atom3 =  None
    #-R: res_in_receptor
    res_in_receptor =  None
    # optional parameters
    ratom1_ratom2_ratom3 = None
    # optional parameters
    verbose = None
    #-g  : rigid output filename
    rigid_filename = None
    #-x  : flexible output filename
    flexres_filename = None

    #'r:l:b:R:B:g:x:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            receptor_filename = a
            #CHECK THAT IT IS ALREADY IN PDBQT FORMAT
            ext = os.path.splitext(receptor_filename)[1]
            if ext!=".pdbqt":
                error = receptor_filename + " not in .pdbqt format"
                raise error
            if verbose: print 'set receptor_filename to ', a
        if o in ('-l', '--l'):
            ligand_filename = a
            ext = os.path.splitext(ligand_filename)[1]
            if ext!=".pdbqt":
                error = ligand_filename + " not in .pdbqt format"
                raise error
            if verbose: print 'set ligand_filename to ', a
        if o in ('-b', '--b'):
            atom1_atom2_atom3 = a
            if verbose: print 'set atom1_atom2_atom3 to ', a
        if o in ('-R', '--R'):
            res_in_receptor = a
            if verbose: print 'set res_in_receptor to ', a
        if o in ('-B', '--B'):
            ratom1_ratom2_ratom3 = a
            if verbose: print 'set ratom1_ratom2_ratom3 to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-g', '--g'):
            rigid_filename = a
            if verbose: print 'set rigid_filename to ', a
        if o in ('-x', '--'):
            flexres_filename = a
            if verbose: print 'set flexres_filename to ', a
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not  receptor_filename:
        print 'prepare_covalent_flexres: receptor filename must be specified!\n'
        usage()
        sys.exit()

    if not  ligand_filename:
        print 'prepare_covalent_flexres: ligand_filename to move must be specified!\n'
        usage()
        sys.exit()

    if not  atom1_atom2_atom3:
        print 'prepare_covalent_flexres: atom1_atom2_atom3 to superimpose must be specified!\n'
        usage()
        sys.exit()

    if not  res_in_receptor:
        print 'prepare_covalent_flexres: res_in_receptor onto which the ligand will superimposed must be specified!\n'
        usage()
        sys.exit()

    if not ratom1_ratom2_ratom3:
        ratom1_ratom2_ratom3 = atom1_atom2_atom3

    #initialize the receptor
    extension = os.path.splitext(receptor_filename)[1]
    if extension!=".pdbqt":
        print 'prepare_covalent_flexres: receptor file must be in .pdbqt format\n'
        usage()
        sys.exit()
    rM = Read(receptor_filename)[0]
    rM.buildBondsByDistance()
    if verbose: print 'read ', receptor_filename

    #initialize the flexible-residue-to-be in the receptor
    resPSet = rM.chains.residues.get(res_in_receptor)
    if not len(resPSet):
        msg= resP.name + " is not a residue in " + receptor_filename
        raise msg
    if len(resPSet)>1:
        msg= resP.name + " specifies more than one residue:" + resP.full_name()
        raise msg
    resP = resPSet[0]
    if verbose: print "receptor residue is %s which has %d atoms" %(resP.full_name(), len(resP.atoms))

    #initialize the ligand
    lM = Read(ligand_filename)[0]
    if verbose: print 'read ', ligand_filename
    lM.buildBondsByDistance()            # try to assign bond order;set to 1 if exception
    #lM.allAtoms.bonds[0].bondOrder = 1   #@@@@@@@@CHECK THIS!!!!@@@@@@@@@
    ah = AtomHybridization()
    ah.assignHybridization(lM.allAtoms)

    #process the ligand:
    #(1)locate specified atom1, atom2 and atom3  SD-CG-CB
    #NB: it is assumed that atom3 has exactly 2 bonds
    at_names = atom1_atom2_atom3.split('_')
    assert len(at_names)==3
    at1Set = lM.allAtoms.get(at_names[0])
    if not len(at1Set):
        msg= at_names[0] + " is not in " + ligand_filename
        raise msg
    at1 = at1Set[0] #SD, closest to ligand, farthest from receptor
    if verbose: print "first ligand atom to superimpose is %s" %(at1.full_name())
    at2Set = lM.allAtoms.get(at_names[1])
    if not len(at2Set):
        msg =  at_names[1] + " is not in " + ligand_filename
        raise msg
    at2 = at2Set[0] #CG
    if verbose: print "second ligand atom to superimpose is %s" %(at2.full_name())
    at3Set = lM.allAtoms.get(at_names[2])
    if not len(at3Set):
        msg =  at_names[2] + " is not in " + ligand_filename
        raise msg
    at3 = at3Set[0] #CB, farthest from ligand, closest to receptor
    if verbose: print "third ligand atom to superimpose is %s" %(at3.full_name())

    #(2)setup shortcut to these two atoms !IN ORDER!
    lig_atoms = AtomSet([at1,at2,at3])
        
    #process the receptor: 
    #(1) locate specified receptor residue
    resSet = rM.chains.residues.get(res_in_receptor)
    if not len(resSet):
        msg = res_in_receptor + " is not in " + receptor_filename
        raise msg
    if len(resSet)!=1:
        msg = res_in_receptor + " selected more than one residue:" + resSet.name
        raise msg
    # SD-CG-CB
    #(2) locate specified ratom1 in that residue
    rat_names = ratom1_ratom2_ratom3.split('_')
    assert len(rat_names)==3
    rat1Set = resP.atoms.get(rat_names[0])
    if not len(rat1Set):
        msg= rat_names[0] + " is not in " + resP.full_name()
        raise msg
    rat1 = rat1Set[0] #SD
    if verbose: print "first receptor atom for superposition is  %s" %(rat1.full_name())
    #(3) locate specified ratom2 in that residue
    rat2Set = resP.atoms.get(rat_names[1])
    if not len(rat2Set):
        msg =  rat_names[1] + " is not in " + resP.full_name()
        raise msg
    rat2 = rat2Set[0] #CG
    if verbose: print "second receptor atom for superposition is  %s" %(rat2.full_name())
    #(4) locate specified ratom3 in that residue
    rat3Set = resP.atoms.get(rat_names[2])
    if not len(rat3Set):
        msg =  rat_names[2] + " is not in " + resP.full_name()
        raise msg
    rat3 = rat3Set[0] #CB, closest to receptor
    if verbose: print "third receptor atom for superposition is  %s" %(rat3.full_name())
    #(5)setup shortcut to these two atoms !IN ORDER!
    ref_atoms = AtomSet([rat1,rat2,rat3])
    if verbose:
        print "for superposition: compute matrix to fit: "
        print "lig_atoms =[at1 = %s, at2 = %s, at3 = %s]" %(at1.full_name(), at2.full_name(), at3.full_name())
        print "onto "
        print "ref_atoms =[rat1 = %s, rat2 = %s, rat3 = %s]" %(rat1.full_name(), rat2.full_name(), rat3.full_name())
    if verbose: print "all input is valid!"
    #(6)locate 'bonds to break' which is the bond from rat3 not to rat2 , here CB: -SD--CG--CB==OtherAtom
    bond_to_break = None
    for b in rat3.bonds:
        otherAt = b.neighborAtom(rat3)
        if otherAt!=rat2:
            bond_to_break = b
            atom_to_attach = otherAt
            if verbose: print "found bond_to_break in receptor => %s" %(b.atom1.name +':'+ b.atom2.name)
            # break the bond
            if verbose: print "receptor atom for ligand attachment in residue %s is %s" %(resP.full_name(),otherAt.full_name())
            otherAt.bonds.remove(b)
            rat3.bonds.remove(b)
            if verbose: print "removed %s-%s bond from %s and from %s"%(b.atom1.name, b.atom2.name, otherAt.full_name(), rat3.full_name())
    if bond_to_break is None:
        print "unable to locate the bond to break in receptor residue %s: expected a bond from %s to an atom other than %s"%(resP.full_name(), rat3.name, rat2.name)
        raise 'invalid input'

    #setup for the superposition
    ref_coords = ref_atoms.coords 
    lig_coords = lig_atoms.coords
    #(1) create the tool
    rigidfitAligner = RigidfitBodyAligner()
    rigidfitAligner.setRefCoords(ref_coords)
    rigidfitAligner.rigidFit(lig_coords)
    #(2) get new coordinates for the moving atoms
    #    here all the ligand atoms move
    new_coords = rigidfitAligner.transformCoords(lM.allAtoms.coords)
    # get the correct shape for the coords (natoms,3)
    #new_coords.shape
    #(25, 4)
    new_coords=new_coords[:,:3] 
    #(3) add the new to ligand
    lM.allAtoms.addConformation(new_coords)
    confNum = len(lM.allAtoms[0]._coords)-1
    #(4) set the coordinates to use to the new ones 
    lM.allAtoms.setConformation(confNum)
    if verbose: print "new coordinates added to ligand"

    #add the ligand atoms other than at1 + at2 to the residue
    #(1) get the atoms to add
    atoms_to_adopt = lM.allAtoms[:]
    #atoms_to_adopt = lM.allAtoms.get(lambda x: not x in lig_atoms)
    if verbose: print "there are %d atoms_to_adopt" %(len(atoms_to_adopt))
    #(3) remove the bonds to the three excluded atoms
    ###attachment_atom.bonds.remove(bond_to_break)
    ###print "removed %s-%s bond from %s"%(bond_to_break.atom1.name, bond_to_break.atom2.name, attachment_atom.full_name())
    ###bond_to_break.neighborAtom(attachment_atom).bonds.remove(bond_to_break)
    if verbose: print 'before removing rat1,rat2 and rat3, residue %s has %d atoms.' %( resP.name, len(resP.atoms))
    #remove rat1, rat2 and rat3
    resP.atoms.remove(rat1)
    if verbose: print 'removed rat1 %s'%rat1.full_name()
    resP.atoms.remove(rat2)
    if verbose: print 'removed rat2 %s'%rat2.full_name()
    resP.atoms.remove(rat3)
    if verbose: print 'removed rat3 %s'%rat3.full_name()
    rM.allAtoms = rM.chains.residues.atoms
    if verbose: print 'after removing rat1,rat2 and rat3 but before adding ligand atoms, residue %s has %d atoms.' %( resP.name, len(resP.atoms))
    if verbose: print 'and the receptor has %d atoms'%len(rM.allAtoms)
    #keep handle to existing atoms for possible clashes
    res_atoms_to_check = resP.atoms[:]
    #(4) add the ligand atoms to the residue
    #for h in atoms_to_adopt:
    for h in lM.allAtoms:
        resP.adopt(h)
        if verbose: print "%s added to %s"%(h.full_name(), resP.full_name())
    if verbose: print 'after adding %d ligand atoms, residue %s has %d atoms.' %( len(lM.allAtoms), resP.name, len(resP.atoms))
    rM.allAtoms = rM.chains.residues.atoms
    if verbose: print 'and the receptor has %d atoms'%len(rM.allAtoms)
    #(5) add the bonds between at3 'CB' in the ligand and the attachment atom in the receptor 
    rM.addBond(at3, atom_to_attach)
    if verbose: 
        print 'built bond between %s from ligand and  atom %s  in residue %s' %( at3.name,atom_to_attach.name, resP.name)
        print "now the bonds of atom %s are:"%(at3.name)
        for b in at3.bonds:
            print " %s -%s"%(b.atom1.full_name(), b.atom2.full_name())

    #(6) do some necessary updates of the protein data structure
    rM.allAtoms = rM.chains.residues.atoms
    num_ats = len(rM.allAtoms)
    rM.allAtoms.number = range(1, num_ats+1)
    # LOOK FOR CLOSE CONTACTS BETWEEN atoms_to_adopt and res_atoms_to_check
    #if find close contacts, iterate rotation around bond
    #in steps of 30degrees using setAngle stuff...
    clashes = detect_clashes(res_atoms_to_check, atoms_to_adopt, verbose=verbose)
    print "@@@ detected %d clashes! @@@"%len(clashes)
    for cl in clashes:
        print "%s - %s = %6.4f"%(cl[0][0].name, cl[0][1].name, cl[1])
    if len(clashes)>0:
        print 'detected %d clash(s)'%(len(clashes))

    # Setup to use AutoDockTools code to write formatted flexres file
    #(1) detect rotatable bonds etc 
    setAutoFlexFields(resP)
    if verbose: print "set autodock fields in atoms in %s" %resP.name
    #to write output
    #(2) setup default output filename if necessary
    #-g  : rigid output filename
    if rigid_filename==None:
        rigid_filename = rM.name + '_rigid.pdbqt'
    if verbose: print "rigid filename is %s" %rigid_filename
    #-x  : flexible output filename
    if flexres_filename==None:
        flexres_filename = rM.name + '_flex.pdbqt'
    if verbose: print "flexible residue filename is %s" %flexres_filename
    #(3) use this class in AutoDockTools for writing formatted outputfiles... 
    frp = AD4FlexibleReceptorPreparation(rM, mode='interactive', 
                rigid_filename=rigid_filename,
                residues=[resP], flexres_filename=flexres_filename)

    frp.write_flex([resP], flexres_filename)
    frp.write_rigid(rM, rigid_filename)
    

# To execute this command type:
# prepare_covalent_flexres.py -l filename -r receptor_filename -l ligand_to_superimpose_filename -b atom1_atom2 -R residue
#specific example:
#python -i prepare_covalent_flexres.py -r 1pwc-protein.pdbqt -l 1pwc-ligand.pdbqt -b OG_CB -R SER62  -g 1pwc_rigid.pdbqt -x 1pwc_flex.pdbqt





