#############################################################################
#
# Author: Ruth HUEY, Stefano FORLI
#
# Copyright: A. Olson TSRI 2010
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/DlgFilters.py,v 1.6 2010/05/07 21:46:28 rhuey Exp $
#
# $Id: DlgFilters.py,v 1.6 2010/05/07 21:46:28 rhuey Exp $
#
#
#
#
#
#
#

"""
Assorted classes of filters for use in selecting subgroups of AutoDock Virtual Screening results...

filter_instance.filter(dlg) returns True/False 

"""
from string import strip
from numpy import oldnumeric as Numeric

class Filter:
    """
        Base class for object to screen dlg results

        input for filter: a dlg instance of the AutoDockTools 'Docking' class
        output: whether the dlg pass the criteria
    """
    def __init__(self, criteria='base_class'): 
        self.criteria = criteria
        self.keyline="BASECLASS"

    
    def get_lines(self, dlg):
        fptr = open(dlg)
        lines = fptr.readlines()
        fptr.close()
        for i in range(len(lines)):
            if lines[i].find(self.keyline)>-1:
                break
        return lines[i:]


    def filter(self, dlg):
        return True



class EnergyFilter(Filter):
    """
        object to evaluate a dlg results based on lowest energy of best docked result

        input: an AutoDockTools.Docking 'dlg' 
        output: whether the best energy of dlg is smaller than energy 
    """
    
    def __init__(self, energy=0):
        self.criteria = "BE"  
        Filter.__init__(self, self.criteria)
        self.energy = energy
        self.keyline= "CLUSTERING HISTOGRAM"

        
    def filter(self, dlg):
        lines = self.get_lines(dlg)
        bestcluster_line = lines[9]
        ll = bestcluster_line.split()
        best_energy = float(ll[2])
        mean_energy = float(ll[6])
        clu_size = int(ll[8])
        return best_energy<self.energy



class ClusterSizeFilter(Filter):
    """
        object to evaluate a dlg results based on size of the largest cluster  

        input: an AutoDockTools.Docking 'dlg' 
        output: whether the largest cluster in the dlg is at least as large as specified cluster_size
    """
    
    def __init__(self, cluster_size=0): 
        self.criteria = "LCsize"
        Filter.__init__(self, self.criteria)
        self.keyline= "CLUSTERING HISTOGRAM"
        self.cluster_size = cluster_size


    def filter(self, dlg):
        lines = self.get_lines(dlg)
        bestenergy_line = lines[9]
        ll = bestenergy_line.split()
        best_energy = float(ll[2])
        largest_clu_size = int(ll[8])
        #check for a larger cluster in lines which follow
        for nl in lines[9:]:
            nl_list = nl.split()
            if len(nl_list)>1 and nl_list[1]=='|':
                this_clu = float(nl_list[8]) 
                if this_clu>largest_clu_size:
                    largest_clu_size = this_clu
            elif nl_list[0]==nl[:-1]:
                break
        return largest_clu_size > self.cluster_size




class ClusterPercentageFilter(Filter):
    """
        object to evaluate a dlg results based on percentage of dlgs in the largest cluster  

        input: an AutoDockTools.Docking 'dlg' 
        output: whether the largest cluster in the dlg contains 
                at least the specified percentage of runs (docked results)
    """
    
    def __init__(self, percentage=0): 
        self.criteria = "LCpercentage"
        Filter.__init__(self, self.criteria)
        self.cluster_percentage = percentage
        self.keyline= "CLUSTERING HISTOGRAM"


    def filter(self, dlg):
        lines = self.get_lines(dlg)
        bestcluster_line = lines[9]
        ll = bestcluster_line.split()
        best_energy = float(ll[2])
        #mean_energy = float(ll[6])
        largest_clu_size = float(ll[8])
        total_clu = 0
        #have to add all the rest to this
        for nl in lines[9:]:
            nl_list = nl.split()
            if len(nl_list)>1 and nl_list[1]=='|':
                this_clu = float(nl_list[8]) 
                total_clu+= this_clu
                if this_clu>largest_clu_size:
                    largest_clu_size = this_clu
            elif nl_list[0]==nl[:-1]:
                break
        percent = 100.*(largest_clu_size/total_clu)
        #print "percent =", percent, " self.cluster_percentage=", self.cluster_percentage
        return percent >= self.cluster_percentage


        
class EnergyClusterSizeFilter(Filter):
    """
        object to evaluate a dlg results based on an energy and a size for the largest cluster
        
        input: an AutoDockTools.Docking 'dlg' 
        output: whether the lowest energy is at least as small as specified value 
                AND 
                whether largest cluster in the dlg is at least as large as specified cluster_size
    """

    def __init__(self, energy=0, cluster_size=0):
        self.criteria = "BE_Csize"  
        Filter.__init__(self, self.criteria)
        self.energyF = EnergyFilter(energy=energy)
        self.cluster_sizeF = ClusterSizeFilter(cluster_size=cluster_size)
        self.cluster_size = cluster_size
        self.keyline= "CLUSTERING HISTOGRAM"


    def filter(self, dlg):
        if not self.energyF.filter(dlg):
            return False
        # so it passed the energy
        if not self.cluster_sizeF.filter(dlg):
            return False
        # so it passed both the energy and cluster_size
        return True


class BestEnergyInLargestCluster(Filter):
    """
        object to screen dlg results based on whether the best energy conformation is in the largest cluster 

        input: an AutoDock dlg 
        output: whether the dlg passes these criteria
    """

    def __init__(self, energy=0):
        self.criteria = "BE_LC"  
        Filter.__init__(self, self.criteria)
        self.keyline= "CLUSTERING HISTOGRAM"


    def filter(self, dlg):
        best_energy = 100
        largest_clu_size = 0
        lines = self.get_lines(dlg)
        bestEnrg = lines[9].split()
        lowest_energy = float(bestEnrg[2])
        le_cluster_size = float(bestEnrg[8])
        #check for a larger cluster size
        ok = True
        for nl in lines[10:]:
            #'_' marks end of cluster output
            if nl[0]=='_': 
                break
            nl_list = nl.split()
            this_clu_size = float(nl_list[8]) 
            if this_clu_size>le_cluster_size:
                largest_clu_size = this_clu_size
                ok = False
                break
        return ok



class EnergyPlusBestEnergyInLargestCluster(Filter):
    """
        object to screen dlg results based on the magnitude of best energy AND 
                                              whether that conf is in the largest cluster 
        input: an AutoDock dlg 
        output: whether the dlg passes these criteria
    """

    def __init__(self, energy=0):
        self.criteria = "BE_LC"  
        Filter.__init__(self, self.criteria)
        self.energyF = EnergyFilter(energy=energy)
        self.be_lcF = BestEnergyInLargestCluster(energy)
        self.keyline= "CLUSTERING HISTOGRAM"


    def filter(self, dlg):
        #check energy first
        if not self.energyF.filter(dlg):
            return False
        #now check whether best energy is in largest cluster
        if not self.be_lcF.filter(dlg):
            return False
        return  True



class LigandEfficiencyFilter(Filter):
    """
        object to screen dlg results based on the 'ligand efficiency' measure which is the 
                            best energy divided by the number of heavy atoms in the ligand

        input: an AutoDock dlg 
        output: whether the dlg passes the criteria
    """

    def __init__(self, ligand_efficiency=0, verbose=False):
        self.criteria = "Ligand_efficiency"  
        Filter.__init__(self, self.criteria)
        self.ligand_efficiency = ligand_efficiency
        self.keyline= "Total number of atoms found"
        self.infoline= "CLUSTERING HISTOGRAM"
        self.h_map_line = ', atom type "HD", grid map index ='
        self.number_of_hydrogens = 0
        self.verbose = verbose

    
    def get_lines(self, dlg):
        fptr = open(dlg)
        lines = fptr.readlines()
        fptr.close()
        atom_numbers = []
        found_index = False
        hmap = -1
        for i in range(len(lines)):
            cur_line = lines[i]
            if not found_index and cur_line.find(self.h_map_line)>-1:
                hmap = int(strip(cur_line.split()[-1]))
                if self.verbose: print "hmap=", hmap
                found_index = 1
            if lines[i].find(self.keyline)>-1:
                break
        for j in range(len(lines[i:])):
            if lines[j].find("Number of atoms with atom type %d ="%(hmap))>-1:
                self.number_of_hydrogens = int(strip(lines[j].split()[-1]))
                if self.verbose: print "found %d hydrogens" %(self.number_of_hydrogens)
        return lines[i:]


    def filter(self, dlg):
        lines = self.get_lines(dlg)
        #only count heavy atoms so need to remove the hydrogens...
        num_lig_atoms = int(lines[0].split()[-2]) - self.number_of_hydrogens
        for i in range(len(lines)):
            if lines[i].find(self.infoline)>-1:
                break
        ll = lines[i+9].split()
        best_energy = float(ll[6])
        lig_eff = float(round(best_energy/num_lig_atoms, 4))
        if self.verbose: print " ligand efficiency = ", lig_eff
        return lig_eff<=self.ligand_efficiency



class InteractionFilter(Filter):       
    """
        class for object to screen dlg results based on presence of specific ligand-receptor interactions
            initialized with a file containing the residue(s) of interest in the receptor

        input: an AutoDock dlg 
        output: whether there are any contacts between a docked pose in the  dlg and the atoms in the 'receptor_file'
    """

    def __init__(self, receptor_file):
        from PyBabel.babelElements import babel_elements
        self.babel_elements = babel_elements
        self.criteria = "Interactions"  
        Filter.__init__(self, self.criteria)
        self.receptor_file = receptor_file
        rptr = open(receptor_file)
        rec_lines = rptr.readlines()
        rptr.close()
        rec_coords =[]
        rec_vdw_rad =[]
        for ll in rec_lines:
            if ll.find("HETATM")>-1 or ll.find("ATOM")>-1:
                rec_coords.append([float(ll[30:38]), float(ll[38:46]),float(ll[46:54])])
                be_key = strip(ll[76:])
                if be_key=='A':
                    be_key = 'C'
                elif be_key=='OA':
                    be_key = 'O'
                elif be_key=='NA':
                    be_key = 'N'
                elif be_key=='SA':
                    be_key = 'S'
                elif be_key=='HD':
                    be_key = 'H'
                rec_vdw_rad.append(babel_elements[be_key]['vdw_rad'])
        # lenC and bigC are about the receptor...
        self.lenC = len(rec_vdw_rad)
        lenC = self.lenC
        #len(rec_vdw_rad)
        #1844
        #len(rec_coords)
        #1844
        checkRadii = Numeric.array(rec_vdw_rad, 'f')
        self.bigRC = Numeric.resize(checkRadii, (lenC, lenC)) #receptor info
        self.bigC = Numeric.resize(rec_coords, (lenC,lenC,3))
        ###bigR = bigRC[:lenK]
        

    def setup_ligand(self, lig_lines):
        lig_coords = []
        first_lig = False
        if not hasattr(self, 'lig_vdw_rad'):
            self.lig_vdw_rad = []
            first_lig = True
        for ll in lig_lines:
            #THESE LINES START WITH 'DOCKED: ' so indices +8
            if ll.find("HETATM")>-1 or ll.find("ATOM")>-1:
                lig_coords.append([float(ll[38:46]), float(ll[46:54]),float(ll[54:62])])
                #lig_coords.append([float(ll[30:38]), float(ll[38:46]),float(ll[46:54])])
                if first_lig: 
                    #indices  +8
                    be_key = strip(ll[84:])
                    if be_key=='A':
                        be_key = 'C'
                    elif be_key=='OA':
                        be_key = 'O'
                    elif be_key=='NA':
                        be_key = 'N'
                    elif be_key=='SA':
                        be_key = 'S'
                    elif be_key=='HD':
                        be_key = 'H'
                    self.lig_vdw_rad.append(self.babel_elements[be_key]['vdw_rad'])
        self.lig_coords = lig_coords
        self.smallM = Numeric.array(lig_coords, 'f')
        #setup variables for ligand arrays
        if first_lig:
            self.lenK=len(lig_coords)
            self.keyRadii = Numeric.array(self.lig_vdw_rad, 'f')
            self.keyRadii.shape = (self.lenK,1)
        self.smallM.shape = (self.lenK,1,3)

        

    def filter(self, dlg):
        fptr = open(dlg)
        dlg_lines = fptr.readlines()
        fptr.close()
        #STEP 1:accumulate lines of various poses
        model_lines = []
        #keep all of them
        all_models = []
        in_model = False
        for ll in dlg_lines:
            if ll.find("DOCKED:")==0:
                #check for a new model
                if ll.find("DOCKED: MODEL")==0:
                    model_lines = []
                in_model = True
                model_lines.append(ll)
            if ll.find("_")==0 and in_model:
                all_models.append(model_lines)
                model_lines = []
                in_model = False
        #initialize this ligand 
        # loop over the models:
        for model_lines in all_models:
            self.setup_ligand(model_lines)
            bigR = self.bigRC[:self.lenK]
            bigM = self.bigC[:self.lenK]
            cutoff = bigR + self.keyRadii
            d = bigM - self.smallM
            dSQ = d*d
            dSQMAT = Numeric.sum(dSQ,2)
            cutoffSQMAT = cutoff*cutoff
            ansMat = Numeric.logical_and(Numeric.less(dSQMAT, cutoffSQMAT),Numeric.not_equal(dSQMAT, 0.))
            rowIndices = Numeric.nonzero(Numeric.sum(ansMat,1))
            num_contacts = 0
            for ind in rowIndices:
                for j in ansMat[ind]: 
                    if j: num_contacts+=1
            if num_contacts > 0:
                break 
        return num_contacts

