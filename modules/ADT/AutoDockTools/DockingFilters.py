#############################################################################
#
# Author: Ruth HUEY, Stefano FORLI
#
# Copyright: A. Olson TSRI 2010
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/DockingFilters.py,v 1.2 2010/04/30 19:53:41 rhuey Exp $
#
# $Id: DockingFilters.py,v 1.2 2010/04/30 19:53:41 rhuey Exp $
#
#
#
#
#
#
#

"""
Assorted classes of filters for use in selecting subgroups of AutoDock Virtual Screening results...

filter_instance.filter(docking) returns True/False 

"""

class Filter:
    """
        Base class for object to screen docking results
        input for filter: a docking instance of the AutoDockTools 'Docking' class
        output: whether the docking pass the criteria
    """
    def __init__(self, criteria='base_class'):
        self.criteria = criteria


    def filter(self, docking):
        return True



class EnergyFilter(Filter):
    """
        object to evaluate a docking results
        based on lowest energy of best docked result
        input: an AutoDockTools.Docking 'docking' 
        output: whether the best energy of docking is smaller than energy 
    """
    
    def __init__(self, energy=0):
        self.criteria = "BE"  
        Filter.__init__(self, self.criteria)
        self.energy = energy


        
    def filter(self, docking, rms=2.0):
        d = docking
        cl = d.clusterer
        clg = cl.clustering_dict[rms]
        val = cl.data[cl.argsort[0]].energy
        return val<self.energy



class ClusterSizeFilter(Filter):
    """
        object to evaluate a docking results
        based on size of the largest cluster  
        input: an AutoDockTools.Docking 'docking' 
        output: whether the largest cluster in the docking is at least as
        large as specified cluster_size
    """
    
    def __init__(self, cluster_size=0): 
        self.criteria = "LCsize"
        Filter.__init__(self, self.criteria)
        self.cluster_size = cluster_size


    def filter(self, docking, rms=2.0):
        d = docking
        clg = d.clusterer.clustering_dict[rms]
        largest_cl = 0
        for cl in clg: 
            #find largest cluster in this docking
            if len(cl)>largest_cl:
                largest_cl = len(cl)
        return largest_cl > self.cluster_size



class ClusterPercentageFilter(Filter):
    """
        object to evaluate a docking results
        based on percentage of dockings in the largest cluster  
        input: an AutoDockTools.Docking 'docking' 
        output: whether the largest cluster in the docking contains 
        at least the specified percentage of dockings
    """
    
    def __init__(self, percentage=0): 
        self.criteria = "LCpercentage"
        Filter.__init__(self, self.criteria)
        self.cluster_percentage = percentage


    def filter(self, docking, rms=2.0):
        d = docking
        clg = d.clusterer.clustering_dict[rms]
        largest_cl = 0
        num_conf = len(d.clusterer.data)
        for cl in clg: 
            #find largest cluster in this docking
            if len(cl)>largest_cl:
                largest_cl = len(cl)
        percent = 100.*(largest_cl/(1.0*num_conf))
        return percent > self.cluster_percentage


        
class EnergyClusterSizeFilter(Filter):
    """
        object to evaluate a docking results based on: 
        an energy and a size for the largest cluster
        input: an AutoDockTools.Docking 'docking' 
        output: whether the lowest energy is at least as small as specified
        value AND whether largest cluster in the docking is at least as
        large as specified cluster_size
    """

    def __init__(self, energy=0, cluster_size=0):
        self.criteria = "BE_Csize"  
        Filter.__init__(self, self.criteria)
        self.energyF = EnergyFilter(energy=energy)
        self.cluster_sizeF = ClusterSizeFilter(cluster_size=cluster_size)
        self.cluster_size = cluster_size


    def filter(self, docking, rms=2.0):
        d = docking
        cl = d.clusterer
        clg = cl.clustering_dict[rms]
        conf = cl.data[cl.argsort[0]]
        e_val = conf.energy
        if not self.energyF.filter(docking,rms):
            return False
        # so it passed the energy
        if not self.cluster_sizeF.filter(docking,rms):
            return False
        # so it passed both the energy and cluster_size
        return True



class EnergyLargestClusterFilter(Filter):
    """
        class for object to screen docking results
        based on the magnitude of best energy and whether that conf is in the largest cluster 
        input: an AutoDock docking 
        output: whether the docking passes these criteria
    """

    def __init__(self, energy=0):
        self.criteria = "BE_LC"  
        Filter.__init__(self, self.criteria)
        self.energyF = EnergyFilter(energy=energy)


    def filter(self, docking, rms=2.0):
        d = docking
        cl = d.clusterer
        clg = cl.clustering_dict[rms]
        conf = cl.data[cl.argsort[0]]
        e_val = conf.energy
        if not self.energyF.filter(docking,rms):
            return False
        #find the largest cluster
        cl_lengths = []
        for cl in clg: cl_lengths.append(len(cl))
        LC_ind = cl_lengths.index(max(cl_lengths))
        largest_cluster = clg[LC_ind]
        #SO is conf in this cluster  or not?
        return conf in largest_cluster



class LigandEfficiencyFilter(Filter):
    """
        class for object to screen docking results
        based on the ligand efficiency measure
        which is the best energy divided by the number 
        of heavy atoms in the ligand
        input: an AutoDock docking 
        output: whether the docking passes the criteria
    """

    def __init__(self, ligand_efficiency=0):
        self.criteria = "Ligand_efficiency"  
        Filter.__init__(self, self.criteria)
        self.ligand_efficiency = ligand_efficiency


    def filter(self, docking):
        d = docking
        cl = d.clusterer
        rms = cl.clustering_dict.keys()[0]
        clg = cl.clustering_dict[rms]
        conf = cl.data[cl.argsort[0]]
        return conf.ligand_efficiency<=self.ligand_efficiency



class InteractionFilter(Filter):       
    """
        @@YET TO BE DONE@@
        class for object to screen docking results
        based on specific interaction ligand-receptor 
        input: an AutoDock docking 
        output: whether the docking passes the criteria
    """

    def __init__(self, interactions=[]):
        self.criteria = "Ligand_efficiency"  
        Filter.__init__(self, self.criteria)
        self.interactions = interactions


    def filter(self, docking):
        d = docking
        cl = d.clusterer
        rms = cl.clustering_dict.keys()[0]
        clg = cl.clustering_dict[rms]
        conf = cl.data[cl.argsort[0]]
        all_passed = False
        # check for presence of all specified interactions
        return all_passed

