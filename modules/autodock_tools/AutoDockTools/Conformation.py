## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Copyright: M. Sanner TSRI 2002
#
# Last modified on Wed May 29 12:28:00 PDT 2002 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Conformation.py,v 1.47 2009/01/09 19:53:49 rhuey Exp $
#

import math
import string
import numpy.oldnumeric as Numeric

from mglutil.math.rmsd import RMSDCalculator
from mglutil.math.statetocoords import StateToCoords
from mglutil.math.transformation import Transformation
from AutoDockTools.cluster import Clusterer
from AutoDockTools.ResultParser import ResultParser



class Conformation:
    """This class can be sent to StateToCoord because attr's match

    """
    def __init__(self, mol, origin, translation, quaternion, torsions,
                 coords = None):
        self.mol = mol
        self.origin = origin
        self.translation = translation
        self.quaternion = quaternion
        self.torsions = torsions
        if coords:
            self.coords = coords
        # the way of the future is (I think) for there to be
        # a state attr and coords attr and a Conformation
        # could be instantiated with either.
        #   do this when stoc is integrated so that stoc.applyState
        #   and deal with it.

        # the cluster dict is keyed by tolerance
        # a conformation knows what cluster(s) it belongs to...
        self.cluster_dict = {}


    def getTorsionOnlyCoords(self):
        """Return your coordinates with no quaterion.

        Don't save these coords, compute them every time.
        """
        self.mol.allAtoms.setConformation(self.mol.stoc.confIndex)
        self.mol.stoc.applyAngList(self.torsions,
            Transformation(trans=self.origin).getMatrix(transpose=1))
        return self.mol.allAtoms.coords


    def getCoords(self):
        """Return your coordinates.

        If the coordinates haven't been computed yet,
            then compute, save, and return them.
        Otherwise, return the previously-computed coordinates.
        """
        #FIX THIS: how could it be None? or []?
        if not hasattr(self, 'coords') or self.coords is None or len(self.coords)==0:
            # then compute the coords
            oldCoords = self.mol.allAtoms.coords[:]
            oldConf = self.mol.allAtoms[0].conformation
            self.mol.allAtoms.setConformation(self.mol.stoc.confIndex)
            self.mol.stoc.applyState(self) # !!! attr's must match !!!
            # and save the coords (this makes a real copy of coords)
            self.coords = Numeric.array(self.mol.allAtoms.coords).tolist()
            self.mol.allAtoms.updateCoords(oldCoords, oldConf)
        return self.coords


    def getRMSD(self, refCoords=None, numCoords=None):
        """Return RMSD of this conformations relative to refCoords.

        If refCoords is not given, the original coordinates for the
        molecule will be used as the reference.
        """
        if not refCoords:
            oldCoords = self.mol.allAtoms.coords[:]
            oldConf = self.mol.allAtoms[0].conformation
            self.mol.allAtoms.setConformation(0)
            refCoords = self.mol.allAtoms.coords[:]
            self.mol.allAtoms.updateCoords(oldCoords, oldConf)
        rmsd_calc = RMSDCalculator(refCoords)            
        if numCoords:
            rmsd = rmsd_calc.computeRMSD(self.getCoords()[:numCoords])
        else:
            rmsd = rmsd_calc.computeRMSD(self.getCoords())
        return rmsd


    def getRMSD_custom(self, refCoords=None):
        """Return the minimum of the regular RMSD and the
        computed RMSD after the coords have been rotated about
        the c2 axis which was aligned with the y-axis.
        """
        normal_RMSD = self.getRMSD(refCoords)

        c2_coords = self.coords[:]
        # I know there's a clever way to use Numeric to do this...
        for c in c2_coords:
            c[0] = -1.0 * c[0]
            c[2] = -1.0 * c[2]

        if not refCoords:
            self.mol.allAtoms.setConformation(0)
            refCoords = self.mol.allAtoms.coords
        rmsd_calc = RMSDCalculator(refCoords)            
        c2_RMSD = rmsd_calc.computeRMSD(self.getCoords())
        return min( normal_RMSD, c2_RMSD)


    def getCoords_subset(self):
        """Return coordinates of current subset, if there is one.

        If the coordinates haven't been computed yet,
            then compute, save, and return them.
        Otherwise, return the previously-computed coordinates.
        """
        if not self.subset:
            return self.getCoords()
        if not hasattr(self, 'subset_coords'):
            oldCoords = self.mol.allAtoms.coords[:]
            oldConf = self.mol.allAtoms[0].conformation
            #by default use the last coords 
            index_to_use = len(self.mol.allAtoms[0]._coords) - 1  
            if hasattr(self.mol, 'stoc'):
                index_to_use = self.mol.stoc.confIndex
            self.mol.allAtoms.updateCoords(self.getCoords(), 
                        index_to_use) # previously: self.mol.stoc.confIndex)
            self.subset_coords = self.subset.coords[:]
            self.mol.allAtoms.updateCoords(oldCoords, oldConf)
        return self.subset_coords


    def getRMSD_subset(self, refCoords=None):
        """Return RMSD of this conformations subset relative to refCoords.

        If refCoords is not given, the original coordinates for the
        subset will be used as the reference.
        """
        if not refCoords:
            refCoords = self.getCoords_subset()
        rmsd_calc = RMSDCalculator(refCoords)            
        return rmsd_calc.computeRMSD(self.getCoords_subset())


    def writeTrj(self, fileptr, istep=1, lastmove='a'):
        """Supply a file_handle and this conformation instance
        will write itself out.
        """
        #FIX THIS: which energy should be written here???
        binding_energy = self.binding_energy
        docking_energy = self.docking_energy
        if not binding_energy:
            binding_energy = 0.000
        if not docking_energy:
            docking_energy = 0.000
        #have to add back the origin to return to AutoDock space
        fileptr.write("state %d %c %f %f %f %f %f %f %f %f %f\n" %
                   (istep,
                    lastmove,
                    binding_energy,
                    docking_energy,
                    self.translation[0] + self.origin[0],
                    self.translation[1] + self.origin[1],
                    self.translation[2] + self.origin[2],
                    self.quaternion[0],
                    self.quaternion[1],
                    self.quaternion[2],
                    self.quaternion[3]))
        for tor in self.torsions:
            fileptr.write("%f\n" % tor)


    def writeRes101(self, fileptr): 
        #data_run_id:
        #1784593 
        #dpf_id:
        #1 
        #ei_version:
        #1.01 
        #ag_version ad_version:
        #3.00 3.05 
        #run_rank: rank of this run in cluster
        #4 
        #run_number: number of this run in cluster
        #6 
        #cluster_rank: rank of this cluster
        #4
        #cluster_size: number of conformations in this cluster
        #1 
        #run_size: number of runs specified in dpf
        #10 
        #rseed1 and rseed2
        #965794405 1006315922 
        #rmsd: RMSD from reference structure
        #26.255011 
        #binding_energy: estimated free energy of binding
        #-8.846132 
        #docking_energy: final docked energy
        #-9.836534 
        #translation x, y, z
        #2.805229 1.033521 1.957878 
        #quaternion unit vector x, y, z, w
        #0.563571 -0.554028 0.612732 -135.133343 
        #number of torsions
        #13 
        #torsion values:
        #-140.26 171.40 -48.43 150.90 -10.93 105.64 -155.73 148.22 174.48 -116.86 -105.52 139.56 -110.30 
        fileptr.write("17 1 1.01 3.00 3.05")
        #nb: conformation.run is the same as run_number
        #nb2:have to add back the origin to return to AutoDock space
        for item in ['run_rank', 'run', 'cluster_rank', 'cluster_size',
                    'run_size', 'rseed1', 'rseed2']:
            if not hasattr(self, item):
                fileptr.write(' 1')
            elif getattr(self, item):
                fileptr.write(" %d" %getattr(self, item))
            else:
                fileptr.write(" 10")
        #write the RMSD for this conformation
        fileptr.write(" %f" %(round(self.getRMSD(),6)))
        for item in [self.binding_energy, self.docking_energy, 
                    self.translation[0] + self.origin[0], 
                    self.translation[1] + self.origin[1], 
                    self.translation[2] + self.origin[2], 
                    self.quaternion[0], self.quaternion[1], 
                    self.quaternion[2], self.quaternion[3]]:
            fileptr.write(" %f" %item)
        fileptr.write(" %d" %len(self.torsions))
        for tor in self.torsions:
            fileptr.write(" %f" % tor)
        #at the end add a newline
        fileptr.write('\n')

    

    def writeRes(self, fileptr):
        #eg
        #9033311 116126 10 1/30/2001 7:27:34 AM 1/30/2001 7:27:34 AM 1.00 3.00 3.05 3 18 3 1 20 980868135 18826794 4.979233 -16.653694 -18.776157 5.005155 -2.431130 15.439371 0.009171 0.617751 0.786320 -170.668379 12 -67.91 -107.72 101.99 -56.92 176.25 21.67 48.10 24.66 -3.65 12.54 -10.66 -3.68 
        #this is for version 1.00
        #1: output_id: %d
        #9033311 
        #2: data_run_id: %d
        #116126 
        #3: dpf_id: %d
        #10 
        #4: creation_dtime: 3part
        #1/30/2001 7:27:34 AM 
        #5: last_update_dtime: 3part
        #1/30/2001 7:27:34 AM 
        #6: ei_version: %f
        #1.00 
        #7: ag_version: %f
        #3.00 
        #8: ad_version: %f
        #3.05 
        #9: run_rank:  %d
        #3 
        #10: run_number:    %d
        #18 
        #11: cluster_rank:  %d
        #3 
        #12: cluster_size:  %d
        #1 
        #13: run_size:  %d
        #20 
        #14: rseed1:    %d
        #980868135 
        #15: rseed2:    %d
        #18826794 
        #16: rmsd:  %f
        #4.979233 
        #a couple of forms of energies: ASK LINDY what these are
        #17: binding_energy: %f
        #-16.653694 
        #18: docking_energy: %f
        #-18.776157 
        #19: translation: 3%f
        #5.005155 -2.431130 15.439371 
        #20: quaternion: 4%f
        #0.009171 0.617751 0.786320 -170.668379 
        #21: number of torsions followed by n values
        #12 -67.91 -107.72 101.99 -56.92 176.25 21.67 48.10 24.66 -3.65 12.54 -10.66 -3.68 
        #write the constant stuff
        fileptr.write("17 18 19 1/23/2001 7:27:34 AM  1/23/2001 7:27:34 AM 1.00 3.00 3.05")
        for item in ['run_rank', 'run', 'cluster_rank', 'cluster_size',
                    'run_size', 'rseed1', 'rseed2']:
            if not hasattr(self, item):
                fileptr.write(' 1')
            elif getattr(self, item):
                fileptr.write(" %d" %getattr(self, item))
            else:
                fileptr.write(" 10")

        #write the RMSD for this conformation
        fileptr.write(" %f" %(round(self.getRMSD(),6)))

        for item in [self.binding_energy, self.docking_energy, 
                    self.translation[0] + self.origin[0], 
                    self.translation[1] + self.origin[1], 
                    self.translation[2] + self.origin[2], 
                    self.quaternion[0], self.quaternion[1], 
                    self.quaternion[2], self.quaternion[3]]:
            fileptr.write(" %f" %item)
        fileptr.write(" %d" %len(self.torsions))
        for tor in self.torsions:
            fileptr.write(" %f" % tor)
        #at the end add a newline
        fileptr.write('\n')


#    def writeRes4(self, fileptr):
#        #eg
#        fileptr.write("17 18 19 1/23/2001 7:27:34 AM  1/23/2001 7:27:34 AM 1.00 3.00 3.05")
#        for item in ['run_rank', 'run', 'cluster_rank', 'cluster_size',
#                    'run_size', 'rseed1', 'rseed2']:
#            if not hasattr(self, item):
#                fileptr.write(' 1')
#            elif getattr(self, item):
#                fileptr.write(" %d" %getattr(self, item))
#            else:
#                fileptr.write(" 10")
#
#        #write the RMSD for this conformation
#        fileptr.write(" %f" %(round(self.getRMSD(),6)))
#
#        for item in [self.binding_energy, self.docking_energy, 
#                    self.translation[0] + self.origin[0], 
#                    self.translation[1] + self.origin[1], 
#                    self.translation[2] + self.origin[2], 
#                    self.quaternion[0], self.quaternion[1], 
#                    self.quaternion[2], self.quaternion[3]]:
#            fileptr.write(" %f" %item)
#        fileptr.write(" %d" %len(self.torsions))
#        for tor in self.torsions:
#            fileptr.write(" %f" % tor)
#        #at the end add a newline
#        fileptr.write('\n')


def dist(c1, c2):
    d = Numeric.array(c2) - Numeric.array(c1)
    ans = math.sqrt(Numeric.sum(d*d))
    return round(ans, 3)

def getAngle(at1, at2, at3 ):
    pt1 = Numeric.array(at1.coords, 'f')
    pt2 = Numeric.array(at2.coords, 'f')
    pt3 = Numeric.array(at3.coords, 'f')
    v1 = Numeric.array(pt1 - pt2)
    v2 = Numeric.array(pt3 - pt2)
    dist1 = math.sqrt(Numeric.sum(v1*v1))
    dist2 = math.sqrt(Numeric.sum(v2*v2))
    sca = Numeric.dot(v1, v2)/(dist1*dist2)
    if sca>1.0:
        sca = 1.0
    elif sca<-1.0:
        sca = -1.0
    ang =  math.acos(sca)*180./math.pi
    return round(ang, 5)

def build_dict(mol):
    d = {}
    bonds = build_bond_dict(mol)
    angles = build_angle_dict(mol)
    d['bonds'] = bonds
    d['angles'] = angles
    return d

def build_bond_dict(mol):
    bonds = {}
    for a in mol.allAtoms:
        for b in a.bonds:
            n1 = b.atom1.number
            n2 = b.atom2.number
            k = (n1,n2)
            if n2<n1:
                k = (n2,n1)
            bonds[k] = dist(b.atom1.coords, b.atom2.coords)
    return bonds


def build_angle_dict(mol):
    angles = {}
    for a1 in mol.allAtoms:
        #print "a1=", a1.name, ':', a1.number
        n1 = a1.number
        #print "n1=",n1
        for b in a1.bonds:
            a2 = b.neighborAtom(a1)
            #print "a2=", a2.name, ':', a2.number
            n2 = a2.number
            #print "n2=",n2
            for b2 in a2.bonds:
                a3 = b2.neighborAtom(a2)
                if a3==a1:
                    continue
                else:
                    n3 = a3.number
                    #print "a3=", a3.name, ':', a3.number
                    k = (a1,a2,a3)
                    kk = (a1.number, a2.number, a3.number)
                    if a1.number>a3.number:
                        kk = (a3.number, a2.number, a1.number)
                    if kk not in angles.keys():
                        new_a = getAngle(a1,a2,a3)
                        angles[kk] = new_a
    return angles
                

class ConformationHandler:
    """This class is bolted onto a Docking instance to manage the conformations which 
    result from an AutoDock experiment, one 'Conformation' per completed run.
    """
    def __init__(self, mol, origin):
        """
        """
        self.mol = mol
        self.origin = origin
        self.conformations = []
        #conformations are 1-based to fit with autodock 1-based rank
        self.current_conf = None
        mol.allAtoms.addConformation(mol.allAtoms.coords)
        self.confIndex = len(mol.allAtoms[0]._coords) - 1
        if hasattr(mol,'torTree'):
            mol.stoc = self.stoc = StateToCoords(mol, origin, self.confIndex)


    def test_conformation(self, conf, confInd, cutoff=0.003, angle_cutoff=15, verbose=False):
        found_different_bond_length = False
        found_different_angle = False
        self.set_conformation(conf)
        test_dict = build_dict(self.mol)
        #test_dict['angles'] = dict
        #test_dict['bonds'] = dict
        for entry in test_dict.keys():
            tD = test_dict[entry]
            rD = self.ref_dict[entry]
            for k,v in tD.items():
                refv = rD[k]
                mag_diff = abs(refv-v)
                if entry=='bonds':
                    if mag_diff>cutoff:
                        found_different_bond_length = True
                        if verbose: print confInd, ": Distance %d%s-%d%s differs: % 6.4f %6.4f \n" \
                            %(k[0]-1,self.mol.allAtoms[k[0]-1].name,k[1]-1, self.mol.allAtoms[k[1]-1].name, refv,v)
                else:
                    if mag_diff>angle_cutoff:
                        found_different_angle = True
                        if verbose: print confInd, ": Angle %d%s-%d%s-%d%s differs: % 6.4f %6.4f \n" \
                            %(k[0]-1,self.mol.allAtoms[k[0]-1].name,k[1]-1, self.mol.allAtoms[k[1]-1].name,k[2]-1, self.mol.allAtoms[k[2]-1].name, refv,v)
        valid = (not found_different_bond_length) and (not found_different_angle)
        if not found_different_bond_length and verbose:
            print confInd, ": no bonds found differing in length more than ", cutoff
        if not found_different_angle and verbose:
            print confInd, ": no angles found differing more than ", angle_cutoff
        return valid
            

    def validate(self, conformations, cutoff=0.003, angle_cutoff=5, verbose=False):
        if len(self.mol.allAtoms[0].bonds)==0:
            self.mol.buildBondsByDistance()
        if not hasattr(self, 'ref_dict'):
            self.ref_dict = build_dict(self.mol)
        valid_conformations = []
        bad_conformations = []
        for ix, conf in enumerate(conformations):
            valid = self.test_conformation(conf,ix,cutoff, angle_cutoff, verbose)
            if valid:
                valid_conformations.append(conf)
            else:
                bad_conformations.append(conf)
        return valid_conformations, bad_conformations
            

    def add(self, clist, keywords=ResultParser.keywords, validate=False, verbose=False, filename=None ):
        """Create/add conformations to the handler.

        clist is a list of dictionaries probably created by
        a subclass of AutoDockTools.ResultParser.
        keywords is the set of keys whose values should become
        Conformation attributes. ResultParser.keywords is the miminal
        and default set.
        """
        #some conformations may not have a translation:
        #in that case want to set translation to (0.,0.,0.)
        #if they do, want to remove it from translation
        #d.get('trn_x', self.origin[0]) returns either trn_x or self.origin[0]
        #thus when you subtract self.origin[0] from what you get you either get
        #corrected translation OR 0. which is what you want... i think
        for d in clist:
            translation = (d.get('trn_x', self.origin[0])- self.origin[0],
                           d.get('trn_y', self.origin[1])- self.origin[1],
                           d.get('trn_z', self.origin[2])- self.origin[2])

            quaternion = (d.get('qtn_nx',0),
                          d.get('qtn_ny',0),
                          d.get('qtn_nz',0),
                          d.get('qtn_ang_deg',0))
            
            newConformation = Conformation(mol = self.mol,
                origin = self.origin,
                # remove origin from AutoDocks reporting of translation
                translation =  translation,
                quaternion =  quaternion,
                torsions = d.get('torsion_values',[]),
                coords = d.get('coords',None))

            newConformation.quaternion0 = (d.get('quaternion_nx',0),
                          d.get('quaternion_ny',0),
                          d.get('quaternion_nz',0),
                          d.get('quaternion_nw',0))

            # add key to newConformation with value from d
            # d is the parser-returned dict
            # missing keys default to None courtesy of d.get
            for key in keywords:
                setattr(newConformation, key, d.get(key, None))
            #NEW: add 'ligand_efficiency' 
            binding_energy = newConformation.binding_energy
            newConformation.ligand_efficiency = 0
            newConformation.filename = filename
            if not hasattr(self.mol, 'lenNonHAtoms'):
                self.mol.lenNonHAtoms = len(self.mol.allAtoms.get(lambda x: x.element!='H'))
            if self.mol.lenNonHAtoms>0:
                newConformation.ligand_efficiency = round(binding_energy/self.mol.lenNonHAtoms, 2)
            # now append to our list
            self.conformations.append(newConformation)
            self.original_conformations = self.conformations

        if validate:
            if verbose: print "pre_validate:len(confs)=", len(self.conformations)
            self.conformations, self.badconformations = self.validate(self.conformations, verbose=verbose)
            if verbose: print "post_validate:len(confs)=", len(self.conformations)


    def set(self, index):
        self.set_conformation(self.conformations[index])


    def set_conformation(self, conf, index=None):
        """Tell the molecule about its new conformation

        """
        # drop the coords into the stoc.confIndex slot of mol._coords
        #self.mol.allAtoms.setConformation(self.confIndex)
        allAtoms = self.mol.allAtoms
        coords = conf.getCoords()
        if not index:
            index = self.confIndex
        allAtoms.updateCoords(coords[:], index)
        #z = self.confIndex
        #for n in range(len(coords)):
        #    allAtoms[n]._coords[z] = coords[n]

        #save possible per atom energy information
        if hasattr(conf, 'estat_energies') and conf.estat_energies and len(conf.estat_energies)==len(allAtoms):
            for n in range(len(allAtoms)):
                allAtoms[n].estat_energy = conf.estat_energies[n]

        if hasattr(conf,'vdw_energies') and conf.vdw_energies and len(conf.vdw_energies)==len(allAtoms):
            for n in range(len(allAtoms)):
                allAtoms[n].vdw_energy = conf.vdw_energies[n]

        if hasattr(conf,'total_energies') and conf.total_energies and len(conf.total_energies)==len(allAtoms):
            for n in range(len(allAtoms)):
                allAtoms[n].total_energy = conf.total_energies[n]

        if hasattr(conf, 'ad4_estat_energies') and conf.ad4_estat_energies and len(conf.ad4_estat_energies)==len(allAtoms):
            for n in range(len(allAtoms)):
                allAtoms[n].ad4_estat_energy = conf.ad4_estat_energies[n]

        if hasattr(conf,'ad4_vdw_energies') and conf.ad4_vdw_energies and len(conf.ad4_vdw_energies)==len(allAtoms):
            #these lump together vdw, hb and dsolv
            for n in range(len(allAtoms)):
                allAtoms[n].ad4_vdw_energy = conf.ad4_vdw_energies[n]

        if hasattr(conf,'ad4_energies') and conf.ad4_energies and len(conf.ad4_energies)==len(allAtoms):
            for n in range(len(allAtoms)):
                allAtoms[n].ad4_energy = conf.ad4_energies[n]

        #conformations are 1-based to fit with autodock 1-based rank
        index = self.conformations.index(conf)
        self.current_conf = conf
        #self.current_conf = index + 1
           

    def writeTrjFile(self, outfptr, run=0, cycle=0, temp=0, istep=0,
                        lastmove=0, eint = 0):
        #first write the header stuff:
        outfptr.write('ntorsions %d\n' % len(self.conformations[0].torsions))
        outfptr.write('run %d\n' % run)
        outfptr.write('cycle %d\n' % cycle)
        outfptr.write('temp %f\n' % temp)
        #what are istep, lastmove and eint????
        istep = 1
        lastmove = 'a'
        for c in self.conformations:
            c.writeTrj(self, outfptr, istep, lastmove)
            istep = istep + 1


    def writeResFile(self, outfptr, istep=0, lastmove=0, eint=0):
        # FIX THIS TO GET id info: istep, lastmove, eint etc
        for c in self.conformations:
            c.writeRes(outfptr, istep, lastmove, eint)
        


class PopulationHandler(ConformationHandler):
    """This class is a specialized ConformationHandler
    designed to handle many individuals in  populations created in a docking.

    """
    def __init__(self, mol, origin):
        """
        """
        ConformationHandler.__init__(self, mol, origin)
        self.individuals = self.conformations
        #individuals conformations are 1-based to fit with autodock 1-based rank
        #self.current_ind = self.conformations
        self.current_pop_ind = 0
        self.all_populations = []
        if len(self.conformations):
            self.all_populations.append(self.conformations)


    def set_current_pop(self, ind):
        assert ind<len(self.all_populations), 'index not in range of populations'
        self.current_pop_ind = ind
        self.conformations = self.all_populations[ind]
        self.individuals = self.conformations
    

    def add(self, clist, keywords=ResultParser.keywords):
        """Create/add conformations to the handler.

        clist is a list of dictionaries probably created by
        a subclass of AutoDockTools.ResultParser.
        keywords is the set of keys whose values should become
        Conformation attributes. ResultParser.keywords is the miminal
        and default set.
        """
        #some conformations may not have a translation:
        #in that case want to set translation to (0.,0.,0.)
        #if they do, want to remove it from translation
        #d.get('trn_x', self.origin[0]) returns either trn_x or self.origin[0]
        #thus when you subtract self.origin[0] from what you get you either get
        #corrected translation OR 0. which is what you want... i think
        confs = []
        for d in clist:
            translation = (d.get('trn_x', self.origin[0])- self.origin[0],
                           d.get('trn_y', self.origin[1])- self.origin[1],
                           d.get('trn_z', self.origin[2])- self.origin[2])

            quaternion = (d.get('qtn_nx',0),
                          d.get('qtn_ny',0),
                          d.get('qtn_nz',0),
                          d.get('qtn_ang_deg',0))

            newConformation = Conformation(mol = self.mol,
                origin = self.origin,
                # remove origin from AutoDocks reporting of translation
                translation =  translation,
                quaternion =  quaternion,
                torsions = d.get('torsion_values',[]),
                coords = d.get('coords',None))

            # add key to newConformation with value from d
            # d is the parser-returned dict
            # missing keys default to None courtesy of d.get
            for key in keywords:
                setattr(newConformation, key, d.get(key, None))

            # now append to our list
            confs.append(newConformation)
        self.all_populations.append(confs)
        self.conformations = self.all_populations[0]
        if len(self.all_populations)==1:
            self.current_pop_ind = 0

    

    

class State:
    """Storage class for the state of any molecule
    """
    def __init__(self, molecule = None, id = 0, origin = (0.,0.,0.),
                 translation = None, quaternion = None,
                 torsions = []):
        """Constructor for class describing molecule State

        self.id = Int              # identifier
        origin = Point()           # center for quaternion motion
        translation = Point()      # translation of ligand center
        quaternion = Quaternion() # orientation of ligand
        torsions = []              # list of torsions in degrees
        """
        self.mol = molecule
        self.id = id
        self.origin = origin
        self.translation = translation
        self.quaternion = quaternion
        self.torsions = torsions
        self.ntorsions = len(torsions)



class AutodockState(State):
    """AdState represents the state of a molecule.
        self.nstep = Int           # number of steps in trajectory cycle
        self.acc_rej_code = ''     # accept/reject code
        e_binding = Float          # intermolecular energy of ligand + macromolecule
                                   # PLUS torsional free energy
        e_total = Float            # energy of ligand + macromolecule
        e_internal = Float         # energy of ligand alone
        e_inter = Float            # Final Intermolecular Energy(from dlg)
        e_intra = Float            # Final Internal Energy of Ligand(from dlg)
        e_tors = Float             # Torsional Free Energy (from dlg)
       The class stores the pdb file name and state
       variables.
    """
    def __init__(self, molecule = None, id = 0, nstep = None, acc_rej_code = None,
                 e_binding = None, e_total = None, e_internal = None, 
                 origin = (0.,0.,0.), translation = None, 
                 quaternion = None, torsions = [], e_inter = None, 
                 e_intra = None, e_tors = None):
        """Constructor for the AutodockState class.
        """
        State.__init__(self, molecule = molecule, id = id, origin = origin,
                translation = translation, quaternion = quaternion, 
                torsions = torsions)
        self.nstep = nstep
        self.acc_rej_code = acc_rej_code
        self.e_binding = e_binding
        self.e_total = e_total
        self.e_internal = e_internal
        self.e_inter = e_inter
        self.e_intra = e_intra
        self.e_tors = e_tors
        assert hasattr(self.mol, 'torTree'), 'molecule has no torTree'
