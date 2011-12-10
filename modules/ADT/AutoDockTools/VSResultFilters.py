#############################################################################
#
# Author: Ruth HUEY, Stefano FORLI
#
# Copyright: A. Olson TSRI 2010
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/VSResultFilters.py,v 1.17 2010/10/05 22:48:37 rhuey Exp $
#
# $Id: VSResultFilters.py,v 1.17 2010/10/05 22:48:37 rhuey Exp $
#
#
#
#
#
#
#

"""
Assorted classes of filters for use in selecting subgroups of AutoDock VSResult files...

filter_instance.filter(VSResult.pdbqt) returns True/False 

10/2010: 
Added a few assorted filters to use in selecting subgroups of AutoDock Vina VSResult files...
"""


class Filter:
    """
        Base class for object to screen AutoDock4 VSResult.pdbqt files

        input for filter: a VSResult.pdbqt instance created by  AutoDockTools/Utilities24/process_VSResults.py 
        output: whether the VSResult passes the criteria
    """
    def __init__(self, criteria='base_class'): 
        self.criteria = criteria
        self.keyline = "BASECLASS"

    
    def get_lines(self, VSResult):
        fptr = open(VSResult)
        lines = fptr.readlines()
        fptr.close()
        #verify that it is a VSResult
        assert lines[0].find('VirtualScreeningResult')>-1
        for i in range(len(lines)):
            if lines[i].find(self.keyline)>-1:
                break
        return lines[i:]


    def filter(self, VSResult):
        return True



class EnergyFilter(Filter):
    """
        object to evaluate an AutoDock4 VSResult  results based on lowest energy of the docked result

        input: an AutoDockTools 'VSResult.pdbqt' 
        output: whether the best energy of VSResult  is smaller than energy 
        ....
        if LE==LC:
        USER  AD> rmsd, LE, clu_size, clu_e_range, dlgfilename, run#, b_curCRDs, b_LE, b_LC
        USER  AD> 0.000,-8.300,92,0.30,faah8930_ZINC00898689_xJ1_xtal_03.dlg,10,1,1,1
        ....

        LE!=LC:
        USER  AD> rmsd, LE, clu_size, clu_e_range, dlgfilename, run#, b_curCRDs, b_LE, b_LC
        USER  AD> 0.000,-5.800,3,1.36,faah8621_ZINC02025973_xJ1_xtal_00.dlg,9,1,1,0
        ....

        if LC:
        USER  AD> rmsd, LE_LC, clu_size, clu_e_range, dlgfilename, run#, b_curCRDs, b_LE, b_LC
        USER  AD> 0.000,-5.800,3,1.36,faah8621_ZINC02025973_xJ1_xtal_00.dlg,9,0,1,0
        USER  AD> 22.435,-5.690,1,0.00,faah8621_ZINC02025973_xJ1_xtal_01.dlg,21,0,0,0
        USER  AD> 15.463,-5.680,2,0.19,faah8621_ZINC02025973_xJ1_xtal_02.dlg,21,0,0,0
        USER  AD> 21.252,-5.590,6,0.69,faah8621_ZINC02025973_xJ1_xtal_00.dlg,14,0,0,0
        USER  AD> 20.951,-5.570,8,0.39,faah8621_ZINC02025973_xJ1_xtal_02.dlg,46,0,0,0
        USER  AD> 21.455,-5.540,2,0.45,faah8621_ZINC02025973_xJ1_xtal_00.dlg,16,0,0,0
        USER  AD> 12.891,-5.480,2,0.42,faah8621_ZINC02025973_xJ1_xtal_01.dlg,14,0,0,0
        USER  AD> 21.837,-5.450,2,0.64,faah8621_ZINC02025973_xJ1_xtal_02.dlg,22,0,0,0
        USER  AD> 5.672,-5.370,2,0.91,faah8621_ZINC02025973_xJ1_xtal_01.dlg,2,0,0,0
        USER  AD> 15.087,-5.340,3,0.93,faah8621_ZINC02025973_xJ1_xtal_02.dlg,1,0,0,0
        USER  AD> 11.694,-5.320,2,0.10,faah8621_ZINC02025973_xJ1_xtal_02.dlg,10,0,0,0
        USER  AD> 5.014,-5.31,9,0.70,faah8621_ZINC02025973_xJ1_xtal_01.dlg,47,1,0,1
        ....


    """
    
    def __init__(self, energy=0):
        self.criteria = "BE"  
        Filter.__init__(self, self.criteria)
        self.keyline = "USER  AD> rmsd,"
        self.energy = energy
        #print "self.energy=", self.energy

        
    def filter(self, VSResult ):
        lines = self.get_lines(VSResult)
        found = 0
        for nl in lines:
            #find current
            #USER  AD> 0.000,-7.470,39,0.60,faah8621_ZINC02026663_xJ1_xtal_00.dlg,30,1,1,1
            ll = nl.strip().split(',')
            #['USER  AD> 0.000', '-7.470', '39', '0.60', 'faah8621_ZINC02026663_xJ1_xtal_00.dlg', '30', '1', '1', '1\n']
            if len(ll)<3:
                break
            if ll[-3]=='1':
                found = 1
                break
        if not found:
            msg =  "improperly formatted file "+ VSResult
            print msg
            return 0
            #raise RuntimeError(msg)
        this_energy = float(ll[1])
        return this_energy<self.energy



class EnergyRangeFilter(Filter):
    """
        object to evaluate an AutoDock4 VSResult  results based on an energy range of best docked result

        input: an AutoDockTools 'VSResult.pdbqt' 
        output: whether the best energy of VSResult is within the specified
        range
        options: inclusive  
            by default inclusive is False so strict comparisons are used. 
            This means that a result passes ONLY if its best_energy is greater than 
            the specified min_energy AND less than the max_energy.
            If inclusive is True, a result passes  if its best_energy is greater than 
                or equal to min_energy AND less than or equal to the specified max_energy
        self.keyline = "USER  AD> binding"
        ....

    """
    
    def __init__(self, min_energy=0, max_energy=0, inclusive=False):
        self.criteria = "BE"  
        Filter.__init__(self, self.criteria)
        #self.keyline = "USER  AD> binding"
        self.keyline = "USER  AD> rmsd,"
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.inclusive = inclusive
        if min_energy>max_energy:
            msg = "minimum energy %f greater than maximum energy %f"%(min_energy, max_energy)
            print msg
            return 0
            #raise RuntimeError(msg)

        
    def filter(self, VSResult):
        found = 0
        lines = self.get_lines(VSResult)
        for nl in lines:
            ll = nl.strip().split(',')
            if len(ll)<3:
                break
            if ll[-3]=='1':
                found = 1
                break
        if not found:
            msg =  "improperly formatted file "+ VSResult
            print msg
            return 0
            #raise RuntimeError(msg)
        best_energy = float(ll[1])
        result = self.min_energy<best_energy<self.max_energy
        if self.inclusive:
            result = self.min_energy<=best_energy<=self.max_energy
        return result




class ClusterSizeFilter(Filter):
    """
        object to evaluate an AutoDock4 VSResult pdbqt file based on size of the cluster  

        input: an AutoDockTools.Docking 'VSResult ' 
        output: whether the largest cluster in the VSResult  is at least as large as specified cluster_size
    """
    
    def __init__(self, cluster_size=0, inclusive=True): 
        self.criteria = "CLsize"
        Filter.__init__(self, self.criteria)
        self.keyline = "USER  AD> rmsd"
        self.cluster_size = cluster_size
        self.inclusive=inclusive

    def filter(self, VSResult ):
        lines = self.get_lines(VSResult)
        found = 0
        for line in lines:
            ll = line.strip().split(',')
            if len(ll)<3:
                break
            if ll[-3]=='1':
                #print "line=", line
                found = 1
                break
        if not found:
            msg =  "improperly formatted file "+ VSResult
            print msg
            return 0
            #raise RuntimeError(msg)
        clu_size = ll[2]
        result = int(clu_size) > self.cluster_size
        if self.inclusive: 
            result = int(clu_size) >= self.cluster_size
        return result



class ClusterPercentageFilter(Filter):
    """
        object to evaluate an AutoDock4 VSResult based on percentage of VSResult  in the largest cluster  

        input: an AutoDockTools.Docking 'VSResult ' 
        output: whether the largest cluster in the VSResult contains 
                at least the specified percentage of runs (docked results)
    """
    
    def __init__(self, percentage=0): 
        self.criteria = "CLpercentage"
        Filter.__init__(self, self.criteria)
        self.keyline = "USER  AD> binding"
        self.cluster_percentage = percentage
        self.verbose = False



    def filter(self, VSResult):
        lines = self.get_lines(VSResult)
        ll = lines[0].split()
        #   0   1      2    3    4   5    6    7
        #USER  AD> binding 2.00 129 runs  31 clusters
        #                        -4   -3  -2   -1
        nruns = int(ll[4])
        found = 0
        for line in lines[2:]:
            #USER  5.014  -5.310    **  9 0.70 *
            ll = line.strip().split(',')
            if len(ll)<3:
                break
            if ll[-3]=='1':
                found = 1
                break
        if not found:
            msg =  "improperly formatted file "+ VSResult
            print msg
            return 0
            #raise RuntimeError(msg)
        this_cl_nruns = int(ll[2])
        if self.verbose: print "this_cl_nruns = ", this_cl_nruns
        percent = 100.*(float(this_cl_nruns)/float(nruns))
        if self.verbose: 
            print "percent =", percent, " self.cluster_percentage=", self.cluster_percentage
        result =  percent >= self.cluster_percentage
        return result


        
class EnergyClusterSizeFilter(Filter):
    """
        object to evaluate an AutoDock4 VSResult results based on an energy and a size for the largest cluster
        
        input: an AutoDockTools.Docking 'VSResult' 
        output: whether the lowest energy is at least as small as specified value 
                AND 
                whether largest cluster in the VSResult is at least as large as specified cluster_size
    """

    def __init__(self, energy=0, cluster_size=0):
        self.criteria = "Energy_ClusterSize"  
        Filter.__init__(self, self.criteria)
        self.energyF = EnergyFilter(energy=energy)
        self.cluster_sizeF = ClusterSizeFilter(cluster_size=cluster_size)
        self.cluster_size = cluster_size


    def filter(self, VSResult):
        if not self.energyF.filter(VSResult):
            return False
        # so it passed the energy
        if not self.cluster_sizeF.filter(VSResult):
            return False
        # so it passed both the energy and cluster_size
        return True



class BestEnergyInLargestCluster(Filter):
    """
        object to screen an AutoDock4 VSResult based on whether the best energy conformation is in the largest cluster 

        input: an AutoDock VSResult 
        output: whether the VSResult passes these criteria
    """

    def __init__(self):
        self.criteria = "BE_LC"  
        Filter.__init__(self, self.criteria)
        self.keyline = "USER  AD> rmsd"

        
    def filter(self, VSResult ):
        lines = self.get_lines(VSResult)
        if lines[0].find("TORSDOF")>-1:
            return False
        lowestenergy_line = lines[1].strip()
        be_is_lc = 0
        if len(lowestenergy_line)>5:
            be_is_lc = lowestenergy_line.find("1,1,1")== len(lowestenergy_line) - 5 
        #be_is_lc = lowestenergy_line.find("1,1,1")!=-1 
        return be_is_lc



class EnergyInLargestCluster(Filter):
    """
        object to screen AutoDock4 VSResult based on the magnitude of best energy in the largest cluster 
        input: an AutoDock VSResult 
        output: whether the VSResult passes these criteria
    """

    def __init__(self, energy=0):
        self.criteria = "Energy_LC"  
        Filter.__init__(self, self.criteria)
        self.energy_criterion = energy
        self.keyline = "USER  AD> rmsd"


    def filter(self, VSResult):
        #check energy first
        lines = self.get_lines(VSResult)
        found = 0
        for line in lines:
            cline = line.strip()
            if cline=='ROOT':
                break
            #if cline.find("1,0,1")>-1 or cline.find("1,1,1")>-1: #??? is there a better way
            #find "USER  AD> 5.014,-5.31,9,0.70,faah8621_ZINC02025973_xJ1_xtal_01.dlg,47,1,0,1"
            #find "USER  AD> 0.000,-7.470,39,0.60,faah8621_ZINC02026663_xJ1_xtal_00.dlg,30,1,1,1"
            #but not "USER  AD> 0.000,-5.930,1,0.00,faah8621_ZINC02058450_xJ1_xtal_00.dlg,1,0,1,0"
            #and not "USER  AD> 0.000,-8.440,10,0.01,faah8621_ZINC02026051_xJ1_xtal_00.dlg,1,1,1,0"
            if cline.find("1,0,1")==len(cline)-5 or cline.find("1,1,1")==len(cline)-5: #??? is there a better way
                found = 1
                break
        if not found:
            #msg =  "improperly formatted file "+ VSResult
            #print msg
            return 0
            #raise RuntimeError(msg)
        #USER  AD> 5.014,-5.31,9,0.70,faah8621_ZINC02025973_xJ1_xtal_01.dlg,47,1,0,1
        ll_list = cline.split(',')
        be_lc = float(ll_list[1])
        #return whether best energy in largest cluster passes...
        return  be_lc<=self.energy_criterion




class EnergyPlusBestEnergyInLargestCluster(Filter):
    """
        object to screen AutoDock4 VSResult based on the magnitude of best energy AND 
                                              whether that conf is in the largest cluster 
        input: an AutoDock VSResult 
        output: whether the VSResult passes these criteria
    """

    def __init__(self, energy=0):
        self.criteria = "BE_LC"  
        Filter.__init__(self, self.criteria)
        self.energyF = EnergyFilter(energy=energy)
        self.be_lcF = BestEnergyInLargestCluster()
        self.keyline = "USER  AD> rmsd"


    def filter(self, VSResult):
        #check energy first
        if not self.energyF.filter(VSResult):
            return False
        #now check whether best energy is in largest cluster
        if not self.be_lcF.filter(VSResult):
            return False
        return  True



class LigandEfficiencyFilter(Filter):
    """
        object to screen AutoDock4 VSResult based on the 'ligand efficiency' measure which is the 
                            best energy divided by the number of heavy atoms in the ligand

        input: an AutoDock VSResult 
        output: whether the VSResult passes the criteria
    """

    def __init__(self, ligand_efficiency=0, verbose=False):
        self.criteria = "Ligand_efficiency"  
        Filter.__init__(self, self.criteria)
        self.ligand_efficiency = ligand_efficiency
        self.keyline = "USER  AD> ligand efficiency"
        self.verbose = verbose


    def filter(self, VSResult):
        lines = self.get_lines(VSResult)
        result = 0
        if lines[0].find('efficiency')>-1:
            ll = lines[0].strip().split()
            lig_eff = float(ll[-1])
            result = lig_eff<=self.ligand_efficiency
        else:
            #@@make this safer
            print VSResult, " missing expected line: 'USER  AD> ligand efficiency float'"
            msg =  "improperly formatted file "+ VSResult
            print msg + ": missing expected line: 'USER  AD> ligand efficiency float'"
            return 0
            #raise RuntimeError(msg)
        return result



class HydrogenBondInteractionFilter(Filter):       
    """
        object to screen AutoDock4 VSResults based on presence of a specific ligand-receptor hydrogen
            bond interaction, initialized with a file containing the residue(s) of interest in the receptor

        input: an AutoDock VSResult 
        output: whether there are any contacts between a docked pose in the  VSResult and the atoms in the 'receptor_file'
    """

    def __init__(self, receptor_str, find_all=False, verbose=False):
        self.criteria = "Interactions"  
        Filter.__init__(self, self.criteria)
        self.keyline = "USER  AD> lig_hb_atoms"
        self.receptor_str = receptor_str
        self.verbose = verbose
        

    def filter(self, VSResult):
        lines = self.get_lines(VSResult)
        #USER lig_hb_atoms : 4
        #USER xJ1_xtal:B:LYS45:NZ,HZ3~ZINC02028532_vs:d:<0>:O2
        #USER xJ1_xtal:B:LYS55:NZ,HZ3~ZINC02028532_vs:d:<0>:O6
        #USER xJ1_xtal:B:MET46:N,HN~ZINC02028532_vs:d:<0>:O4
        #USER xJ1_xtal:B:LYS55:NZ,HZ3~ZINC02028532_vs:d:<0>:O7
        num_hb_atoms = int(lines[0].split()[-1])
        found_interaction = False
        for i in range(num_hb_atoms):
            if lines[i+1].find(self.receptor_str)>-1:
                found_interaction = True
                break  #
        return found_interaction



class CloseContactInteractionFilter(Filter):       
    """
        object to screen AutoDock4 VSResult based on presence of specific ligand-receptor interactions
            initialized with a file containing the residue(s) of interest in the receptor

        input: an AutoDock VSResult 
        output: whether there are any contacts between a docked pose in the  VSResult and the atoms in the 'receptor_file'
    """

    def __init__(self, receptor_str, find_all=False, verbose=False):
        self.criteria = "Interactions"  
        Filter.__init__(self, self.criteria)
        self.keyline = "USER  AD> macro_close_ats:"
        self.receptor_str = receptor_str
        self.verbose = verbose
        

    def filter(self, VSResult):
        lines = self.get_lines(VSResult)
        #USER macro_close_ats: 19
        num_contacts = int(lines[0].split()[-1])
        found_interaction = False
        for i in range(num_contacts):
            if lines[i+1].find(self.receptor_str)>-1:
                found_interaction = True
                break  #
        return found_interaction



class VinaFilter:
    """
        Base class for object to screen VinaVSResult.pdbqt files
        input for filter: a VinaVSResult.pdbqt file created with AutoDockTools/Utilities24/process_VinaResult.py 
        output: whether the VinaVSResult passes the criteria
    """
    def __init__(self, criteria='base_class'): 
        self.criteria = criteria
        self.keyline = "BASECLASS"

    
    def get_lines(self, VSResult):
        fptr = open(VSResult)
        lines = fptr.readlines()
        fptr.close()
        #verify that it is a VSResult
        assert lines[0].find('AutoDockVina VirtualScreeningResult')>-1
        for i in range(len(lines)):
            if lines[i].find(self.keyline)>-1:
                break
        return lines[i:]


    def filter(self, VSResult):
        return True



class VinaEnergyFilter(VinaFilter):
    """
        object to evaluate a VinaVSResult  results based on energy 

        input: an AutoDockTools 'VinaVSResult.pdbqt' 
        output: whether the best energy of VinaVSResult  is smaller than energy 
        ....
        USER  AD>  2 of 9 MODELS
        USER  AD>  ligand efficiency  -0.2217
        USER  AD>  -5.10, 26.865, 27.587
        ....

    """
    
    def __init__(self, energy=0):
        self.criteria = "Energy"  
        VinaFilter.__init__(self, self.criteria)
        self.keyline = "USER  AD>  ligand efficiency"
        self.energy = energy

        
    def filter(self, VinaVSResult ):
        lines = self.get_lines(VinaVSResult)
        #lines[0]
        #USER  AD>  ligand efficiency  -0.2217
        #lines[1]
        #USER  AD>  -5.10, 26.865, 27.587
        info_line = lines[1] 
        ll = info_line.split()
        if len(ll)==5 and ll[3][-1]==',':
            found = 1
        else:
            msg =  "improperly formatted file "+ VSResult
            print msg
            return 0
            #raise RuntimeError(msg)
        this_energy = float(ll[2][:-1])
        #print " this_energy=", this_energy
        return this_energy<self.energy



class VinaEnergyRangeFilter(VinaFilter):
    """
        object to evaluate a VinaVSResult  results based on an energy range of best docked result
        input: an AutoDockTools 'VinaVSResult.pdbqt' 
        output: whether the energy of VinaVSResult is within the specified range
        options: inclusive  
            by default inclusive is False so strict comparisons are used. 
            This means that a result passes ONLY if its best_energy is greater than 
            the specified min_energy AND less than the max_energy.
            If inclusive is True, a result passes  if its best_energy is greater than 
                or equal to min_energy AND less than or equal to the specified max_energy
        self.keyline = "USER  AD>  ligand efficiency"
        ....

    """
    
    def __init__(self, min_energy=0, max_energy=0, inclusive=False):
        if min_energy>max_energy:
            msg = "minimum energy %f greater than maximum energy %f"%(min_energy, max_energy)
            print msg
            return 0
        self.criteria = "Energy"  
        VinaFilter.__init__(self, self.criteria)
        self.min_filter = VinaEnergyFilter(min_energy)
        self.max_filter = VinaEnergyFilter(max_energy)
        self.inclusive = inclusive


    def filter(self, VinaVSResult ):
        min_ok = self.min_filter.filter(VinaVSResult)
        max_ok = self.max_filter.filter(VinaVSResult)
        #print "returning ", max_ok and not min_ok
        return max_ok and not min_ok


class VinaLigandEfficiencyFilter(VinaFilter):
    """
       object to screen VinaVSResult results based on the 'ligand efficiency' measure which is the 
                           best energy divided by the number of heavy atoms in the ligand

       input: an AutoDock VinaVSResult.pdbqt 
       output: whether the VinaVSResult passes the criteria
    """
    def __init__(self, ligand_efficiency=0, verbose=False):
        self.criteria = "Ligand_efficiency"  
        VinaFilter.__init__(self, self.criteria)
        self.ligand_efficiency = ligand_efficiency
        self.keyline = "USER  AD>  ligand efficiency"
        self.verbose = verbose


    def filter(self, VinaVSResult ):
        lines = self.get_lines(VinaVSResult)
        #lines[0]
        #USER  AD>  ligand efficiency  -0.2217
        #lines[1]
        #USER  AD>  -5.10, 26.865, 27.587
        info_line = lines[0] 
        ll = info_line.split()
        if len(ll)==5:
            found = 1
        else:
            msg =  "improperly formatted file "+ VinaVSResult
            print msg
            return 0
            #raise RuntimeError(msg)
        this_efficiency = float(ll[4].strip())
        #print " this_efficiency =", this_efficiency
        return this_efficiency<self.ligand_efficiency


###NO VinaHydrogenBondInteractionFilter because spherical model used for hydrogen positions
###class HydrogenBondInteractionFilter(Filter):       
###    """
###        class for object to screen VSResult results based on presence of a specific ligand-receptor hydrogen
###            bond interaction, initialized with a file containing the residue(s) of interest in the receptor

###        input: an AutoDock VSResult 
###        output: whether there are any contacts between a docked pose in the  VSResult and the atoms in the 'receptor_file'
###    """


class VinaCloseContactInteractionFilter(VinaFilter):       
    """
        class for object to screen VinaVSResult results based on presence of specific ligand-receptor interactions
            initialized with a file containing the residue(s) of interest in the receptor

        input: an AutoDock VinaVSResult 
        output: whether there are any contacts between a docked pose in the  VSResult and the atoms in the 'receptor_file'
    """

    def __init__(self, receptor_str, find_all=False, verbose=False):
        self.criteria = "Interactions"  
        VinaFilter.__init__(self, self.criteria)
        self.keyline = "USER  AD> macro_close_ats:"
        self.receptor_str = receptor_str
        self.verbose = verbose
        

    def filter(self, VinaVSResult):
        lines = self.get_lines(VinaVSResult)
        #USER macro_close_ats: 12
        num_contacts = int(lines[0].split()[-1])
        found_interaction = False
        for i in range(num_contacts):
            if lines[i+1].find(self.receptor_str)>-1:
                found_interaction = True
                break  #only need to find first...
        return found_interaction


