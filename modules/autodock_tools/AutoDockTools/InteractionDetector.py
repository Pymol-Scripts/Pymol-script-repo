#############################################################################
#
# Author: Ruth HUEY, Stefano FORLI
#
# Copyright: A. Olson TSRI 2010
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/InteractionDetector.py,v 1.17 2010/10/06 18:07:19 rhuey Exp $ 
#
# $Id: InteractionDetector.py,v 1.17 2010/10/06 18:07:19 rhuey Exp $
#
#
#
#
#
#
#

"""
Assorted classes for use in detecting interactions of various subgroups in AutoDock Virtual Screening results...

interactionDetector_instance.screen(pdbqt_result) returns interactions_found as True/False. 
If interactions_exist, a new pdbqt_result file including information describing the interactions is written.  
Types of interactions to be detected include:
    -contacts between receptor atoms and ligand atoms closer than sum of vdw radii
    -hydrogen-bond interactions between receptor atoms and ligand atoms of appropriate types
    -pi-pi interactions between receptor atoms and ligand atoms of appropriate types
    -cation-pi interactions between receptor atoms and ligand atoms of appropriate types
    -t_shaped interactions between receptor atoms and ligand atoms of appropriate types

The filter is initialized with a file containing the receptor residues whose interactions with ligands are of interest.
For example, in the case of hsg1, a file could potentially contain only A:ARG8 and B:ARG8 atoms.

"""
from string import strip
from numpy import oldnumeric as Numeric
from MolKit import Read
from MolKit.molecule import MoleculeSet, Atom, AtomSet, HydrogenBond, BondSet
from MolKit.protein import Residue, ResidueSet
from MolKit.distanceSelector import CloserThanVDWSelector, DistanceSelector
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder
from PyBabel.cycle import RingFinder
from PyBabel.aromatic import Aromatic
from MolKit.bondSelector import AromaticCycleBondSelector
import numpy.oldnumeric as Numeric
from mglutil.math import crossProduct


class InteractionDetector:
    """
        Base class for object to detect interactions between a receptor and potential ligands... such as virtual screen results
        
        initialized with the receptor file containing entire molecule or residues of interest.
        processLigand method takes as input a pdbqt file containing docked coordinates and returns a string 
                        composed of  hbondStr + macrocloseContactStr + ligcloseContactStr
        If interections are found, a new pdbqt containing interaction description is also output: 
    """
    def __init__(self, receptor_file, percentCutoff=1., detect_pi=False, dist_cutoff=6, verbose=False, 
                        distanceSelector=None, hydrogen_bond_builder=None, distanceSelectorWithCutoff=None,
                        aromatic_cycle_bond_selector=None): 
        self.receptor_file = receptor_file
        receptor = Read(receptor_file)
        assert(len(receptor)==1)
        assert isinstance(receptor, MoleculeSet)
        self.macro = receptor[0]
        self.macro_atoms = self.macro.allAtoms
        self.macro.buildBondsByDistance()
        self.verbose = verbose
        #??useful??
        self.percentCutoff = percentCutoff
        self.detect_pi = detect_pi
        self.distanceSelector = distanceSelector 
        if self.distanceSelector is None:
            self.distanceSelector = CloserThanVDWSelector(return_dist=0)
        self.hydrogen_bond_builder = hydrogen_bond_builder 
        if self.hydrogen_bond_builder is None:
            self.hydrogen_bond_builder = HydrogenBondBuilder()
        self.distanceSelectorWithCutoff = distanceSelectorWithCutoff 
        if self.distanceSelectorWithCutoff is None:
            self.distanceSelectorWithCutoff = DistanceSelector()
        self.dist_cutoff=float(dist_cutoff)
        self.results = d = {}
        self.report_list =['lig_hb_atoms','lig_close_atoms']
        self.aromatic_cycle_bond_selector = aromatic_cycle_bond_selector
        if self.aromatic_cycle_bond_selector is None:
            self.aromatic_cycle_bond_selector = AromaticCycleBondSelector()
        if detect_pi:
            self.report_list.extend(['pi_cation','pi_pi', 'cation_pi', 't_shaped'])
        if self.verbose: print "self.report_list=", self.report_list


    def processLigand(self, ligand_file, percentCutoff=None, output=1, outputfilename=None, comment="USER  AD> ", buildHB=1, remove_modelStr=False):
        #if outputfilename is not None, have to (1) update all the strings or (2) rename ligand
        self.results = d = {}
        ligand = Read(ligand_file)
        assert isinstance(ligand, MoleculeSet)
        assert(len(ligand)==1)
        ligand = ligand[0]
        first = ligand.name.find('_model')
        last = ligand.name.rfind('_model')
        if first!=last:
            ligand.name = ligand.name[:last]
        if remove_modelStr:
            curName = ligand.name
            modelIndex=curName.find("_model")
            if modelIndex>-1:
                curName = curName[:modelIndex]
            #setup outputfilestem
            ligand.name = curName
        ligand.buildBondsByDistance()
        if not percentCutoff:
            percentCutoff = self.percentCutoff
        # first detect sets of atoms forming hydrogen bonds
        # hbondStr,hblist and has_hbonds added to process vina results which do not include valid hydrogen bonds
        hbondStr = ""
        hblist = [""]
        has_hbonds = False
        if buildHB:
            has_hbonds = True
            hbondStr = self.buildHydrogenBonds(ligand, comment=comment)           #
            if self.verbose: print "hbondStr=", hbondStr
            hblist = hbondStr.split('\n') # hbond info
            if hblist[0]=='lig_hb_atoms : 0':
                has_hbonds = False
        # next detect sets of atoms in close contact not forming hydrogen bonds
        macrocloseContactStr, ligcloseContactStr = self.buildCloseContactAtoms(percentCutoff, ligand, comment=comment)              #
        self.piResults = ""
        if self.detect_pi: 
            self.detectPiInteractions()
            if self.results['pi_cation']:
                self.piResults = self.get_pi_cation_result(print_ctr=1, comment=comment)
                if self.verbose: 
                    print "found pi_cation in ", ligand_file
                    print "set self.piResults to ", self.piResults
            if self.results['pi_pi']:
                self.piResults += self.get_pi_pi(print_ctr=1, comment=comment)
                if self.verbose: 
                    print "set self.piResults to ", self.piResults
        macro_cclist = macrocloseContactStr.split(';')
        lig_cclist = ligcloseContactStr.split(';')
        if has_hbonds or len(macro_cclist):
            fptr = open(ligand_file)
            lines = fptr.readlines()
            fptr.close()
            if outputfilename is not None:
                optr = open(outputfilename, 'w')
            else:
                optr = open(ligand_file, 'w')
            for new_l in hblist:
                if not(len(new_l)): #don't output empty lines
                    continue
                if new_l.find(comment)!=0:
                    new_l = comment + new_l
                if new_l != hblist[-1]:
                    optr.write(new_l+"\n")
                else:
                    optr.write(new_l) # no newline after the last one
            for new_l in macro_cclist:
                if not(len(new_l)): #don't output empty lines
                    continue
                if new_l.find(comment)!=0:
                    new_l = comment + new_l
                if new_l != macro_cclist[-1]:
                    optr.write(new_l + "\n")
                else:
                    optr.write(new_l)
            for new_l in lig_cclist:
                if not(len(new_l)): #don't output empty lines
                    continue
                if new_l.find(comment)!=0:
                    new_l = comment + new_l
                if new_l != lig_cclist[-1]:
                    optr.write(new_l + "\n")
                else:
                    optr.write(new_l)
            if len(self.piResults): ##???
                for s in self.piResults.split("\n"):
                    optr.write(s+"\n")
            for l in lines:
                optr.write(l)
            optr.close()
        return hbondStr + macrocloseContactStr + ligcloseContactStr
            

    def getResultStr(self, verbose=False):
        #only write important ones: hbonds, vdw_contacts, ?pi_pi,pi_cation etc?
        ss = ""
        for k in self.report_list:
            v = self.results[k]
            #key:atom1.full_name();
            if len(v):
                if verbose: print "gRS: report for ", k
                if verbose: print  "USER:  " + k + ":" + str(len(v)) 
                if verbose: print  "     start ss= ", ss
                ss +=  k + ":" + str(len(v))  + "\n"
                if k=='lig_hb_atoms':
                    if verbose: print "@@ getResultStr lig_hb_atoms"
                    hbstr = ""
                    for lig_at in v:
                        for hb in lig_at.hbonds:
                            if hasattr(hb, 'hAt'):
                                hbstr += 'USER %s-%s~%s\n'%(hb.donAt.full_name(), hb.hAt.full_name(),hb.accAt.full_name())
                            else:
                                hbstr += 'USER %s~%s\n'%(hb.donAt.full_name(), hb.accAt.full_name())
                    if verbose: print "hbstr="
                    #add it to ss here
                    ss += hbstr
                    if verbose: print "with hbstr: ss=", ss
                #if self.verbose:
                elif k == 'pi_cation':
                    for res in v:
                        #v=[(<Atom instance> HSG1:A:ARG8:CZ, [<Atom instance> IND: : INI 20:C1])]
                        #res=(<Atom instance> HSG1:A:ARG8:CZ, [<Atom instance> IND: : INI 20:C1])
                        #res[0]=<Atom instance> HSG1:A:ARG8:CZ
                        #res[1]=  [<Atom instance> IND: : INI 20:C1]
                        #res[1][0]=  <Atom instance> IND: : INI 20:C1
                        #
                        ss += "USER %s~~%s\n"%(res[0].full_name(), res[1][0].full_name())
                else:
                    for w in v:
                        if verbose: print "@@ getResultStr w=", w
                        try:
                            ss += 'USER ' + w.full_name()+'\n;'
                        except:
                            if verbose: print "except on ", w
                            ss += "USER " + str(w)
                ss += "\n" 
                if verbose: print  "     end ss= ", ss
        return ss


    def buildHydrogenBonds(self, ligand, comment="USER AD> "):
        h_pairDict = self.hydrogen_bond_builder.build(ligand.allAtoms, self.macro_atoms)
        self.h_pairDict = h_pairDict
        #keys should be from lig, values from macro 
        #sometimes are not...@@check this@@
        h_results = {}
        for k, v in h_pairDict.items():
            h_results[k] = 1
            for at in v:
                h_results[at] = 1
        all_hb_ats = AtomSet(h_results.keys())  #all
        d = self.results
        macro_hb_ats = d['macro_hb_atoms'] = all_hb_ats.get(lambda x: x.top==self.macro)
        self.macro_hb_ats = macro_hb_ats
        # process lig
        lig_hb_ats = d['lig_hb_atoms'] = all_hb_ats.get(lambda x: x in ligand.allAtoms)
        self.lig_hb_ats = lig_hb_ats
        outS = comment + "lig_hb_atoms : %d\n"%(len(lig_hb_ats))
        for p in self.lig_hb_ats:  #intD.results['lig_hb_atoms']:
            for hb in p.hbonds:
                if hasattr(hb, 'used'): continue
                if hb.hAt is not None:
                    outS += comment + "%s,%s~%s\n"%(hb.donAt.full_name(), hb.hAt.name, hb.accAt.full_name())
                else:
                    outS += comment + "%s~%s\n"%(hb.donAt.full_name(), hb.accAt.full_name())
                hb.used = 1
                #hsg1V:B:ARG8:NH2,HH22~clean: : INI 20:N5
                #clean: : INI 20:O4,H3~hsg1V:B:ASP29:OD2
                #clean: : INI 20:O4,H3~hsg1V:B:ASP29:OD2
                #clean: : INI 20:N4,H3~hsg1V:B:GLY27:O
                #clean: : INI 20:O2,H2~hsg1V:B:ASP25:OD1
        #macroHStr = self.macro.allAtoms.get(lambda x: hasattr(x, 'hbonds') and len(x.hbonds)).full_name()
        #ligHStr = ligand.allAtoms.get(lambda x: hasattr(x, 'hbonds') and len(x.hbonds)).full_name()
        #return  macroHStr + '==' + ligHStr
        if self.verbose: 
            print  "buildHB returning:"
            print outS
        return outS



    def buildCloseContactAtoms(self, percentCutoff, ligand, comment="USER AD> "):
        pairDict = self.distanceSelector.select(ligand.allAtoms,
                        self.macro_atoms, percentCutoff=percentCutoff)
        self.pairDict = pairDict
        #reset here
        lig_close_ats = AtomSet()
        macro_close_ats = AtomSet()
        cdict = {}
        for k,v in pairDict.items():
            if len(v):
                cdict[k] = 1
            for at in v:
                if at not in macro_close_ats:
                    cdict[at] = 1
        closeAtoms = AtomSet(cdict.keys())
        lig_close_ats = closeAtoms.get(lambda x: x.top==ligand).uniq()
        #ligClAtStr = lig_close_ats.full_name()
        ligClAtStr = comment + "lig_close_ats: %d\n" %( len(lig_close_ats))
        if len(lig_close_ats):
            ligClAtStr += comment + "%s\n" %( lig_close_ats.full_name())
        macro_close_ats = closeAtoms.get(lambda x: x in self.macro_atoms).uniq()
        macroClAtStr = comment + "macro_close_ats: %d\n" %( len(macro_close_ats))
        if len(macro_close_ats):
            macroClAtStr += comment + "%s\n" %( macro_close_ats.full_name())
        #macroClAtStr = "macro_close_ats: " + len(macro_close_ats)+"\n" +macro_close_ats.full_name()
        rdict = self.results
        rdict['lig_close_atoms'] = lig_close_ats
        rdict['macro_close_atoms'] = macro_close_ats
        if self.verbose: print "macroClAtStr=", macroClAtStr
        if self.verbose: print "ligClAtStr=", ligClAtStr
        if self.verbose: print "returning "+ macroClAtStr + '==' + ligClAtStr
        return  macroClAtStr , ligClAtStr


    def getCations(self, atoms):
        #select atoms in ARG and LYS residues
        arg_cations = atoms.get(lambda x: (x.parent.type=='ARG' and \
                                x.name in ['CZ']))
        lys_cations = atoms.get(lambda x: (x.parent.type=='LYS' and \
                                x.name in ['NZ', 'HZ1', 'HZ2', 'HZ3']))
        #select any positively-charged metal ions... cannot include CA here
        metal_cations = atoms.get(lambda x: x.name in ['Mn','MN', 'Mg',\
                                'MG', 'FE', 'Fe', 'Zn', 'ZN'])
        ca_cations = atoms.get(lambda x: x.name in ['CA', 'Ca'] and x.parent.type=='CA')
        cations = AtomSet() 
        #cations.extend(arg_cations)
        for a in arg_cations:
            cations.append(a)
        #cations.extend(lys_cations)
        for a in lys_cations:
            cations.append(a)
        #cations.extend(metal_cations)
        for a in metal_cations:
            cations.append(a)
        #cations.extend(ca_cations)
        for a in ca_cations:
            cations.append(a)
        return cations


    def detectPiInteractions(self, tolerance=0.95, debug=False, use_all_cycles=False):
        if debug: print "in detectPiInteractions"
        self.results['pi_pi'] = []        #stacked rings...?
        self.results['t_shaped'] = []     #one ring perpendicular to the other
        self.results['cation_pi'] = []    #
        self.results['pi_cation'] = []    #
        self.results['macro_cations'] = []#
        self.results['lig_cations'] = []  #
        #at this point have self.results
        if not len(self.results['lig_close_atoms']):
            return
        lig_atoms = self.results['lig_close_atoms'].parent.uniq().atoms
        macro_res = self.results['macro_close_atoms'].parent.uniq()
        if not len(macro_res):
            return
        macro_atoms = macro_res.atoms
        l_rf = RingFinder()
        #Ligand
        l_rf.findRings2(lig_atoms, lig_atoms.bonds[0])
        #rf.rings is list of dictionaries, one per ring, with keys 'bonds' and 'atoms'
        if debug: print "LIG: len(l_rf.rings)=", len(l_rf.rings)
        if not len(l_rf.rings):
            if debug: print "no lig rings found by l_rf!"
            return
        acbs = self.aromatic_cycle_bond_selector
        #acbs = AromaticCycleBondSelector()
        lig_rings = []
        for r in l_rf.rings:
            ring_bnds = r['bonds']
            if use_all_cycles: 
                lig_rings.append(ring_bnds)
            else:
                arom_bnds = acbs.select(ring_bnds)
                if len(arom_bnds)>4:
                    lig_rings.append(arom_bnds)
        if debug: print "LIG: len(lig_arom_rings)=", len(lig_rings)
        self.results['lig_rings'] = lig_rings
        self.results['lig_ring_atoms'] = AtomSet()
        #only check for pi-cation if lig_rings exist
        if len(lig_rings):
            macro_cations = self.results['macro_cations'] = self.getCations(macro_atoms)
            macro_cations = macro_cations.get(lambda x: x.element!='H')
            lig_ring_atoms = AtomSet()
            u = {}
            for r in lig_rings:
                for a in BondSet(r).getAtoms():
                    u[a] = 1
            if len(u): 
                lig_ring_atoms = AtomSet(u.keys())
                lig_ring_atoms.sort()
                self.results['lig_ring_atoms'] = lig_ring_atoms
            if len(macro_cations):
                if debug: print "check distances from lig_rings to macro_cations here"
                #macro cations->lig rings
                pairDict2 = self.distanceSelector.select(lig_ring_atoms,macro_cations)
                z = {}
                for key,v in pairDict2.items():
                    val = v.tolist()[0]
                    if val in macro_cations:
                        z[val] = [key]
                if len(z):
                    self.results['pi_cation'] = (z.items())
                else:
                    self.results['pi_cation'] = []
        #check the distance between the rings and the macro_cations
        self.results['lig_cations'] = self.getCations(lig_atoms)
        lig_cations = self.results['lig_cations'] 
        #remove hydrogens
        lig_cations = lig_cations.get(lambda x: x.element!='H')
        #Macromolecule
        m_rf = RingFinder()
        m_rf.findRings2(macro_res.atoms, macro_res.atoms.bonds[0])
        #rf.rings is list of dictionaries, one per ring, with keys 'bonds' and 'atoms'
        if debug: print "MACRO: len(m_rf.rings)=", len(m_rf.rings)
        if not len(m_rf.rings):
            if debug: print "no macro rings found by m_rf!"
            return
        macro_rings = []
        for r in m_rf.rings:
            ring_bnds = r['bonds']
            if use_all_cycles: 
                macro_rings.append(ring_bnds)
            else:
                arom_bnds = acbs.select(ring_bnds)
                if len(arom_bnds)>4:
                    macro_rings.append(arom_bnds)
        if debug: print "len(macro_arom_rings)=", len(macro_rings)
        self.results['macro_rings'] = macro_rings
        self.results['macro_ring_atoms'] = AtomSet()
        #only check for pi-cation if macro_rings exist
        if len(macro_rings):
            macro_ring_atoms = AtomSet()
            u = {}
            for r in macro_rings:
                for a in BondSet(r).getAtoms(): #new method of bondSets
                    u[a] = 1
            if len(u):
                macro_ring_atoms = AtomSet(u.keys())
                macro_ring_atoms.sort()
                self.results['macro_ring_atoms'] = macro_ring_atoms
            if len(lig_cations):
                if debug: print "check distances from macro_rings to lig_cations here"
                pairDict3 = self.distanceSelector.select(macro_ring_atoms,lig_cations)
                z = {}
                for x in pairDict3.items():
                    #lig cations->macro rings
                    z.setdefault(x[1].tolist()[0], []).append(x[0])
                if len(z):
                    self.results['cation_pi'] = (z.items())
                else:
                    self.results['cation_pi'] = []
                #macro_pi_atoms = AtomSet(pairDict3.keys())
                #l_cations = AtomSet()
                #for v in pairDict3.values():
                #    for x in v:
                #        l_cations.append(x)
                #self.results['cation_pi'] = pairDict3.items()
                #self.results['cation_pi'] = (l_cations, macro_pi_atoms)
        #check for intermol distance <6 Angstrom (J.ComputChem 29:275-279, 2009)
        #compare each lig_ring vs each macro_ring
        for lig_ring_bnds in lig_rings:
            lig_atoms = acbs.getAtoms(lig_ring_bnds)
            lig_atoms.sort()
            if debug: print "len(lig_atoms)=", len(lig_atoms)
            #---------------------------------
            # compute the normal to lig ring
            #---------------------------------
            a1 = Numeric.array(lig_atoms[0].coords)
            a2 = Numeric.array(lig_atoms[2].coords)
            a3 = Numeric.array(lig_atoms[4].coords)
            if debug: print "a1,a2, a3=", a1.tolist(), a2.tolist(), a3.tolist()
            for macro_ring_bnds in macro_rings:
                macro_atoms = acbs.getAtoms(macro_ring_bnds)
                macro_atoms.sort()
                if debug: print "len(macro_atoms)=", len(macro_atoms)
                pD_dist = self.distanceSelectorWithCutoff.select(macro_ring_atoms, lig_atoms, cutoff=self.dist_cutoff)
                if not len(pD_dist[0]):
                    if debug: 
                        print "skipping ligand ring ", lig_rings.index(lig_ring_bnds), " vs ",
                        print "macro ring", macro_rings.index(macro_ring_bnds)
                    continue
                #---------------------------------
                # compute the normal to macro ring
                #---------------------------------
                b1 = Numeric.array(macro_atoms[0].coords)
                b2 = Numeric.array(macro_atoms[2].coords)
                b3 = Numeric.array(macro_atoms[4].coords)
                if debug: print "b1,b2, b3=", b1.tolist(), b2.tolist(), b3.tolist()
                # check for stacking 
                a2_1 = a2-a1
                a3_1 = a3-a1
                b2_1 = b2-b1
                b3_1 = b3-b1
                if debug: print "a2_1 = ", a2-a1
                if debug: print "a3_1 = ", a3-a1
                if debug: print "b2_1 = ", b2-b1
                if debug: print "b3_1 = ", b3-b1
                n1 = crossProduct(a3_1,a2_1) #to get the normal for the first ring
                n2 = crossProduct(b3_1,b2_1) #to get the normal for the second ring
                if debug: print "n1=", n1
                if debug: print "n2=", n2
                n1 = Numeric.array(n1)
                n2 = Numeric.array(n2)
                n1_dot_n2 = Numeric.dot(n1,n2)
                if debug: print "n1_dot_n2", Numeric.dot(n1,n2)
                if abs(n1_dot_n2) >= 1*tolerance: 
                    if debug: print "The rings are stacked vertically" 
                    new_result = (acbs.getAtoms(lig_ring_bnds), acbs.getAtoms(macro_ring_bnds))
                    self.results['pi_pi'].append(new_result)
                if abs(n1_dot_n2) <= 0.01*tolerance: 
                    if debug: print "The rings are stacked perpendicularly" 
                    new_result = (acbs.getAtoms(lig_ring_bnds), acbs.getAtoms(macro_ring_bnds))
                    self.results['t_shaped'].append(new_result)


    def get_pi_pi(self, print_ctr=1, comment="USER AD> "):
        #one ring parallel to another 
        #@@ UNTESTED! ... need test data
        if not len(self.results.get( 'pi_pi')):
            return ""
#        self.results['pi_pi'] = []        #stacked rings...?
        res = self.results['pi_pi']
        ss = comment + "pi_pi: %d\n"%(len(self.results['pi_pi']))
        ####
        #3l1m_lig_vs:A:HEM150:C4A,C3D,C1A,C2D,C4D,C2A,NA,ND,C3A,C1D : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:C3B,NB,C4A,C1A,C4B,NA,C2A,C2B,C3A,C1B : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:C3B,C1B,C3C,NB,C4C,NC,C2B,C1C,C4B,C2C : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:NC,C4C,C1C,C3C,C2C : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:C1A,NA,C2A,C3A,C4A : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:NB,C2B,C3B,C1B,C4B : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:NC,C4C,C1C,C3C,C2C : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #3l1m_lig_vs:A:HEM150:C4D,ND,C1D,C3D,C2D : 3l1m_rec:A:PHE65:CD2,CD1,CG,CZ,CE1,CE2
        #
        # build string for  self.results['pi_pi'] 
        for res in self.results['pi_pi']:
            ss += comment + "%s~~%s\n"%(res[0].full_name(), res[1].full_name())
            ###ss += "USER %s~~%s\n"%(res[0].full_name(), res[1].full_name())
            #ss += "USER %s~~%s\n"%(res[0].full_name(), res[1][0].full_name())
        #print "returning ss=" , ss
        return ss


    def get_t_shaped(self, print_ctr=1, comment="USER AD> "):
        #one ring perpendicular to the other
        #@@ UNTESTED! ... need test data
        if not len(self.results.get( 't_shaped')):
            return ""
        #print "in gpcr: self.results[t_shaped]=", self.results['t_shaped']
        ss = comment + "t_shaped: %d\n"%(len(self.results['t_shaped']))
        # build string for  self.results['t_shaped'] 
        for res in self.results['t_shaped']:
            #for example:
            #????
            #????
            #????
            #
            ss += comment + "%s~~%s\n"%(res[0].full_name(), res[1][0].full_name())
        return ss


    def get_cation_pi(self, print_ctr=1, comment="USER AD> "):
        #cation near aromatic ring 
        #@@ UNTESTED! ... need test data
        if not len(self.results.get( 'cation_pi')):
            return ""
        ss = comment + "cation_pi: %d\n"%(len(self.results['cation_pi']))
        # build string for  self.results['cation_pi'] 
        for res in self.results['cation_pi']:
            #for example:
            #????
            #????
            #????
            #
            ss += comment + "%s~~%s\n"%(res[0].full_name(), res[1][0].full_name())
        return ss


    def get_pi_cation_result(self, print_ctr=1, comment="USER AD> "):
        #aromatic ring near cation 
        #print "in gpcr: self.results[pi_cation]=", self.results['pi_cation']
        if not len(self.results.get( 'pi_cation')):
            return ""
        ss = comment + "pi_cation: %d\n"%(len(self.results['pi_cation']))
        # build string for self.results['pi_cation'] 
        for res in self.results['pi_cation']:
            #for example:
            #v=[(<Atom instance> HSG1:A:ARG8:CZ, [<Atom instance> IND: : INI 20:C1])]
            #res=(<Atom instance> HSG1:A:ARG8:CZ, [<Atom instance> IND: : INI 20:C1])
            #res[0]=<Atom instance> HSG1:A:ARG8:CZ
            #res[1]=  [<Atom instance> IND: : INI 20:C1]
            #res[1][0]=  <Atom instance> IND: : INI 20:C1
            #
            #attempt to get names for all atoms in the whole ring:
            res1_str = res[1][0].full_name()
            for rr in self.results['lig_rings']:
                ats = rr.getAtoms()
                if res[1][0] in ats:
                    res1_str = ats.full_name()
                    break
            ss += comment + "%s~~%s\n"%(res[0].full_name(), res1_str)
            #ss += "USER %s~~%s\n"%(res[0].full_name(), res[1][0].full_name())
        return ss



    def print_ligand_residue_contacts(self, print_ctr=1):
        ctr = 1
        #pairDict is atom-based 
        #for residue-based report:
        # need to build lists of unique parents of keys 
        # and lists of unique parents of corresponding values
        res_d = {}
        for at in self.pairDict.keys():
            if at.parent not in res_d.keys():
                res_d[at.parent] = {}
            for close_at in self.pairDict[at]:
                res_d[at.parent][close_at.parent] = 1
        #print it out
        for lig_res in res_d.keys():
            if print_ctr:
                print ctr, lig_res.parent.name+':'+ lig_res.name + '->',
            else:
                print lig_res.parent.name+':'+ lig_res.name + '->',
            for macro_res in res_d[lig_res]:
                print macro_res.parent.name + ':' + macro_res.name + ',',
            print
            ctr += 1
        return res_d


    def print_macro_residue_contacts(self, print_ctr=1):
        ctr = 1
        #pairDict is atom-based 
        #for residue-based report:
        # need to build lists of unique parents of keys 
        # and lists of unique parents of corresponding values
        res_d = {}
        for at_key, at_list in self.pairDict.items():
            for at in at_list:
                if at.parent not in res_d.keys():
                    res_d[at.parent] = {}
                res_d[at.parent][at_key.parent] = 1
        #print it out
        for macro_res in res_d.keys():
            if print_ctr:
                print ctr, macro_res.parent.name+':'+ macro_res.name + '->',
            else:
                print macro_res.parent.name+':'+ macro_res.name + '->',
            for lig_res in res_d[macro_res]:
                print lig_res.parent.name + ':' + lig_res.name + ',',
            print
            ctr += 1
        return res_d

    def print_report(self, keylist=[]):
        if not len(keylist):
            keylist = [
                'lig_close_atoms',
                'lig_hb_atoms',
                'lig_hbas',
                'macro_close_atoms',
                'macro_hb_atoms',
                    ]
        d = self.results
        if self.verbose:
            for k in keylist:
                print k, ':', len(d[k]), '-', d[k].__class__


    def print_hb_residue(self, print_ctr=1):
        ctr = 1
        #pairDict is atom-based 
        res_d = {}
        for at in self.h_pairDict.keys():
            if at.parent not in res_d.keys():
                res_d[at.parent] = {}
            for close_at in self.h_pairDict[at]:
                res_d[at.parent][close_at.parent] = 1
        # print it out
        # Hbond instance are define with donAt,hAt ~ accAtt  (tilde)
        for don_res in res_d.keys():
            if print_ctr:
                print ctr, don_res.top.name+':'+don_res.parent.name+':'+ don_res.name + '->',
            else:
                print don_res.top.name+':'+don_res.parent.name+':'+ don_res.name + '->',
            for acc_res in res_d[don_res]:
                print acc_res.top.name+':'+acc_res.parent.name + ':' + acc_res.name + ',',
            print
            ctr += 1
        return res_d

