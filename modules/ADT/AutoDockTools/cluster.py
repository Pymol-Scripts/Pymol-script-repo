## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Mon Mar 15 17:47:37 PST 2004 by lindy
#
# $Id: cluster.py,v 1.41 2010/06/21 15:47:17 rhuey Exp $
#
"""
Some definitions:
    cluster:    a list of conformations all within tolerance RMSD.
    clustering: a list of clusters from clustering a set of conformations
                at a given tolerance RMSD.
    Clusterer:  instances of this class perform the clustering operation,
                maintain the distance_matrix, and keep the dictionary of
                clusterings (keyed by tolerance).
"""

import string
import sys
import numpy.oldnumeric as Numeric
import UserList

class Cluster(UserList.UserList):
    """A cluster is a list of conformations.

    The 'rank' of a conformation is its order in the list. 
    """
    def __init__(self, seed, info=None):
        UserList.UserList.__init__(self, [seed])
        self.seed = seed
        if info:
            self.build(info)


    def build(self, info):
        #info is a dictionary 
        #confs is a ordered list of Conformations
        self.confs = info['confs']
        self.do_stats(info['min_energy'],
                      info['max_energy'],
                       info['average_energy'])
        self.rank = info['rank']
   

    def do_stats(self, min_energy=None, max_energy=None, avg_energy=None):
        # set me
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.avg_energy = avg_energy



class Clustering(UserList.UserList):
    """A clustering is list of clusters resulting from a clustering operation.
    """
    def __init__(self):
        UserList.UserList.__init__(self)
        self.tolerance = None # set by the Clusterer


    def do_stats(self):
        """Compute some statistics about this clustering.
        """
        self.avg_size = Numeric.sum(map(len, self))/float(len(self))
        for c in self:
            c.do_stats()


        
class Clusterer:
    """This class knows how to cluster a list of conformations
    the Autodock way. That is sort conformations by energy, visit
    every conformation in order adding to the cluster if within the
    tolerance RMSD. Seed a new cluster other wise.

    self.data: the list of conformations to cluster
    self.clustering: a float-keyed dictionary. Keys are clustering
    tolerances (in angstroms).
    
    The list of Cluster instances is self.clusters
    """
    def __init__(self, data, sort='binding'):
        """data is a list of Conformations.
        """
        self.data = data
        # extract the binding (default) or docking energy
        if sort == 'docking':
            energy_list = [d.docking_energy for d in data]
        elif sort == 'intermolecular':
            energy_list = [d.intermol_energy for d in data]
        elif sort == 'flexres':
            if hasattr(data[0], 'subset'):
                ind = len(data[0].getCoords()) - len(data[0].subset) 
            else:
                print "no subsets have been set up for conformations!"
                return "ERROR"
            for d in data:
                d.flexres_energy = Numeric.add.reduce(d.total_energies[ind:])
                d.flexres_index = ind
            energy_list = [d.flexres_energy for d in data]
        else:
            energy_list = [d.binding_energy for d in data]

        # sort the conformations by energy
        self.argsort = Numeric.argsort(energy_list)
        self.energy_used = sort

        # save the pair-wise distances for reuse by get_distance
        self.dist_matrix = Numeric.zeros([len(data), len(data)]) - 1.0
        # set the customizable get_distance method to default
        self.set_get_distance(self._get_distance_default)
        
        self.clustering_dict = {}


    def _get_distance_flexres(self, a, b):
        """return RMSD between two conformations, a and b.
        """
        ax = self.data.index(a)
        bx = self.data.index(b)
        assert a.flexres_index == b.flexres_index
        ind = a.flexres_index

        if self.dist_matrix[ax][bx] >= 0.0:
            # return previously saved distance
            return self.dist_matrix[ax][bx]
        else:
            # compute, save, and return distance
            dist = a.getRMSD_subset(b.getCoords()[b.flexres_index:])
            self.dist_matrix[ax][bx] = self.dist_matrix[bx][ax] = dist
            return dist


    def _get_distance_custom(self, a, b):
        """return RMSD between two conformations, a and b.
        """
        ax = self.data.index(a)
        bx = self.data.index(b)

        if self.dist_matrix[ax][bx] >= 0.0:
            # return previously saved distance
            return self.dist_matrix[ax][bx]
        else:
            # compute, save, and return distance
            dist = a.getRMSD_custom(b.getCoords())
            self.dist_matrix[ax][bx] = self.dist_matrix[bx][ax] = dist
            return dist


    def _get_distance_default(self, a, b):
        """return RMSD between two conformations, a and b.
        """
        ax = self.data.index(a)
        bx = self.data.index(b)

        if self.dist_matrix[ax][bx] >= 0.0:
            # return previously saved distance
            return self.dist_matrix[ax][bx]
        else:
            # compute, save, and return distance
            dist = a.getRMSD(b.getCoords())
            self.dist_matrix[ax][bx] = self.dist_matrix[bx][ax] = dist
            return dist


    def _get_distance_subset(self, a, b):
        """return RMSD between subsets of atoms in two conformations, a and b.
        """
        ax = self.data.index(a)
        bx = self.data.index(b)

        if self.dist_matrix[ax][bx] >= 0.0:
            # return previously saved distance
            return self.dist_matrix[ax][bx]
        else:
            ind = len(a.getCoords()) - len(a.subset)
            # compute, save, and return distance
            dist = a.getRMSD_subset(b.getCoords()[ind:])
            self.dist_matrix[ax][bx] = self.dist_matrix[bx][ax] = dist
            return dist


    def set_get_distance(self, f):
        self.get_distance = f


    def get_distance(self, a, b):
        return self.get_distance


    def make_clustering(self, tolerance, ref=None):
        """
        tolerance is the RMSD tolerance in angstroms
        ref is the reference Conformation
        """
        this_clustering = Clustering()
        # seed the first cluster with supplied ref
        # or the lowest energy conformation in the data (as ref)
        if not ref:
            ref = self.data[int(self.argsort[0])]
        this_clustering.append( Cluster(ref))

        # go through all conformations
        for i in range(len(self.argsort)):
            conf = self.data[int(self.argsort[i])]
            if conf == ref: continue  # seen it before so skip it

            # try to find a home for this conf in all the this_clustering
            clustered = 0
            for cluster in this_clustering:
                # get dist from cluster seed to conformaton
                dist = self.get_distance(cluster[0], conf)

                if dist <= tolerance:
                    # then add to current cluster and stop looking
                    cluster.append(conf)
                    clustered = 1
                    break
            if not clustered:
                # then start new cluster
                this_clustering.append( Cluster(conf))

        # save the clustering in the dictionary
        self.clustering_dict[tolerance] = this_clustering
        this_clustering.tolerance = tolerance

        # tell each conformation its cluster and rank (for this tolerance)
        for cx, cluster in enumerate(this_clustering):
            for rank, conf in enumerate(cluster):
                conf.cluster_dict[tolerance] = (cx, rank)


    def get_max_dist(self, clust1, clust2):
        max_dist = -1.
        for conf_i in clust1:
            for conf_j in clust2:
                dist_ij = self.get_distance(conf_i, conf_j)
                if dist_ij>max_dist: max_dist = dist_ij
        return max_dist


    def get_min_dist(self, clust1, clust2):
        min_dist = 1000000000.
        for conf_i in clust1:
            for conf_j in clust2:
                dist_ij = self.get_distance(conf_i, conf_j)
                if dist_ij<min_dist: min_dist = dist_ij
        return min_dist


    def make_hierarchical_clustering(self, linkage='single', precision=1, debug=False):
        """ uses HierarchicalClustering of python-cluster as explained in 
        http://www.elet.polimi.it/upload/matteucc/Clustering/tutorial_html/hierarchical.html
        Available linkage choices are: 'single','complete', 'average' or 'uclus'"""
        try:
            from python_cluster.cluster import HierarchicalClustering
        except:
            print "HierarchicalClustering from python-cluster.cluster not found"
            return
        assert linkage in ['single', 'complete', 'average', 'uclus']
        cl = HierarchicalClustering(self.data, self.get_distance, linkage=linkage)
        #make sure all the distances have been calculated
        num_confs = len(self.data)
        for i in range(num_confs):
            for j in range(num_confs):
                self.get_distance(self.data[i], self.data[j])
        all_dists = self.dist_matrix.ravel().tolist()
        all_dists.sort()
        #remove the zeros from the diagonal
        dists = all_dists[num_confs:] 
        dists.sort()
        rmsds = []
        clustD = self.clustering_dict
        duplicates = []
        #precision parameter was added to try to permit
        #the user to control how many levels were calculated...
        # by default, round rmsd values to 1 decimal place
        # eg 0.1, 0.2 etc...
        #@@ need to test with data for which this might matter
        # with ind.dlg, it doesn't seem to make any difference...
        for i in dists:
            v = round(i, precision)
            if v not in rmsds:
                rmsds.append(v)
        for i in range(len(rmsds)):            
            level = rmsds[i]
            if debug: print "getlevel level=", level
            clustD[level] = cl.getlevel(level)
            if debug: print "built->len(clustD[",level,"]=", len(clustD[level])
            if len(clustD[level])==1:
                if debug: print "all in a single cluster"
  

    def show_clustering(self):
        dk = self.clustering_dict.keys()
        dk.sort()
        for val in dk:
            print val, ':'
            for clust in self.clustering_dict[val]:
                for conf in clust:
                    print conf.run,
                print
            print        
                    
    def write_summary(self, filename=None):
        """
        """
        if filename:
            file_ptr = open( filename, 'w')
        else:
            file_ptr = sys.stdout

        # make sure the keys happen in the same order
        sorted_keys = self.clustering_dict.keys()[:]
        sorted_keys.sort()
        # write tolerances
        for tolerance in sorted_keys:
            clustering = self.clustering_dict[tolerance]
            clustering.do_stats()
            file_ptr.write("%s  %3d  %3.2f" % (str(tolerance),
                                           len(clustering),
                                           clustering.avg_size))
            file_ptr.write("\n")


    def write(self, filename=None):
        """Write a set of clusterings to a file.

The clusterings are written to the filename given in
the following format:
1st line: space separated floats describing tolerances.
the rest: one line foreach conformation in order of
          increasing binding_energy. each line consiststs
          of space separated integers, two per clustering,
          the first gives the cluster index, the second
          gives the rank within the cluster.

For example:
1.0 2.0 3.0
0 0 0 0 0 0
0 1 0 2 0 0

New parameter for docking/binding energy denotes which was
used to create clusterings:
For example:
1.0 2.0 3.0 docking  -- specifies that DOCKING energy was used
1.0 2.0 3.0 binding  -- specifies that BINDING energy was used
1.0 2.0 3.0          -- if not specified then BINDING energy assumed.

Mon Mar 15 16:09:51 PST 2004
New fields reporting the rmsd values from the reference strucuture
(or lowest energy cluster seed if no reference structure given) and
rmsd values to each of the cluster seeds. Finally the energy of the
conformations is given. This is the docking energy or binding energy
as specified in the header. (This energy is for convenience as it is
just a copy of the energy from the AutoDock output.

For example:
1.5 2.0 binding
  0   0   0   0  0.000  0.000  0.000 -14.645
  1   0   1   0  2.449  0.000  0.000 -14.636
  0   1   0   1  1.281  1.281  1.281 -14.424
  1   1   1   1  2.548  1.014  1.014 -14.210
                   |      |      |
                   |      |      |-rmsd from cluster seed @ 2.0 tolerance
                   |      |-rmsd from cluster seed @ 1.0 tolerance
                   |-rmsd from overall reference structure      

The first line says there were three clusterings at 1.5 and 2.0 Angstrom
rmsd tolerances using the binding energy as the sort. The second conformation
started a new cluster because it was 2.449 A-rmsd from the the reference.
        """
        if filename:
            file_ptr = open( filename, 'w')
        else:
            file_ptr = sys.stdout

        # make sure the keys happen in the same order
        sorted_keys = self.clustering_dict.keys()[:]
        sorted_keys.sort()
        # write tolerances
        for key in sorted_keys:
            file_ptr.write("%6s  " % str(key))
        # note which energy was used to cluster
        file_ptr.write(self.energy_used)
        file_ptr.write("\n")

        # write one line per conformation
        for i in range(len(self.argsort)):
            conf = self.data[int(self.argsort[i])]
            # write the pair for each tolerance key
            for key in sorted_keys:
                file_ptr.write("%3d %3d " % conf.cluster_dict[key])
                
            # append the rmsd from the overall reference conf (or seed)
            file_ptr.write("%6.3f " % self.get_distance(
                self.clustering_dict[key][0][0],
                conf))
            # append the rmsd from the clustering seeds
            for key in sorted_keys:
                file_ptr.write("%6.3f " % self.get_distance(
                    self.clustering_dict[key][conf.cluster_dict[key][0]][0],
                    conf))
            # and finally the appropriate energy
            if self.energy_used=='energy':
                file_ptr.write("%7.3f " %
                           getattr(conf, "binding_energy" ))
            else:
                file_ptr.write("%7.3f " %
                           getattr(conf, self.energy_used + "_energy" ))
            # end the line
            file_ptr.write("\n")
        if file_ptr != sys.stdout:
            file_ptr.close()


    def getInfoStr(self, comment='USER  AD>', ind=-1, rms=-1, ncl_to_write=1, include_max=False, report_all=False, include_dlgfilename_run=False):
        """Write a set of clusterings to a string instead of file 

A string summarizing the clusterings is returned with the format as described in
write method above:

Default:
binding\n 0   0   0   0  0.000  0.000  0.000 -14.645\n 1   0   1   0  2.449  0.000  0.000 -14.636\n 0   1   0   1  1.281  1.281  1.281 -14.424\n 1   1   1   1  2.548  1.014  1.014 -14.210\n

If comment is '#':
#binding\n#  0   0   0   0  0.000  0.000  0.000 -14.645\n#  1   0   1   0  2.449  0.000  0.000 -14.636\n#  0   1   0   1  1.281  1.281  1.281 -14.424\n#  1   1   1   1  2.548  1.014  1.014 -14.210\n

        """
        first = comment
        #if comment: first=comment
        clu_str = "%s"%first
        ind = int(ind)
        #corresponding filenames+ run numbers
        dlg_run_str = ""

        # make sure the keys happen in the same order
        sorted_keys = self.clustering_dict.keys()[:]
        sorted_keys.sort()
        if rms==-1:
            rms = sorted_keys[0]
            print "reporting rms %f clustering:"%rms
        if rms not in sorted_keys:
            print "no clustering exists at rms %f"%rms
            return "ERROR"
        # add tolerances and energy used on 1 line
        if report_all and len(sorted_keys)>1:
            for key in sorted_keys[1:]:
                clu_str += "%6s  " % str(key)
        # note which energy was used to cluster
        clu_str += "%6s\n" %str(self.energy_used)
        # add one line per cluster
        out_cl_num = len(self.clustering_dict[rms])
        n_omitted_cl = 0
        n_omitted_confs =0
        e_range_omitted_cl = 0
        nconf_omitted_cl = 0
        if ncl_to_write<out_cl_num:
            out_cl_num = ncl_to_write
            #record information about clusters not in infoStr
            #number of omitted clusters
            len_dict = len(self.clustering_dict[rms])
            n_omitted_cl = len_dict - ncl_to_write # output
            #print "set n_omitted_cl to ", n_omitted_cl
            #number of conformations
            #for cl in self.clustering_dict[rms][ncl_to_write:]:
            #print "len_dict-ncl_to_write=", len_dict - ncl_to_write
            for j in range(ncl_to_write, len_dict):
                cl = self.clustering_dict[rms][j]
                n_omitted_confs += len(cl)                # output
                #print "adding ", len(cl), " now n_omitted_confs=", n_omitted_confs
            #first and last omitted cluster energy
            first_omitted_e = self.clustering_dict[rms][ncl_to_write][0].binding_energy
            #??should last energy be best of last cluster or worst of last cluster??
            last_omitted_e = self.clustering_dict[rms][-1][-1].binding_energy
        for i in range(out_cl_num):
            #conf is lowest energy conformation in 'i-th' cluster
            conf = self.clustering_dict[rms][i][0]
            # append the rmsd from the overall reference conf (or seed)
            clu_str+="%s%.3f," %(first, self.get_distance(
                self.clustering_dict[rms][0][0],
                conf))
            # and finally the appropriate energy
            if self.energy_used=='energy':
                clu_str+="%.3f" % getattr(conf, "binding_energy" )
            else:
                clu_str+="%.3f" % getattr(conf, self.energy_used + "_energy" )
            # check whether ind corresponds to this conf
            #if ind==int(self.argsort[i]):
            if include_dlgfilename_run:
                dlg_run_str += "%s,%d" %(conf.filename, conf.run)
            #if self.data[ind]==conf:
            #    #clu_str += " * %d %d"%(ind, i)
            #    dlg_run_str += "1"
            #    #clu_str += "1,"
            #else:
            #    dlg_run_str += "0"
            ## end line for this conformation
            clu_str+="\n"
            dlg_run_str += "\n"
        if n_omitted_cl:
            clu_str="USER omitted %d clusters [%d confs]: be range %6.4f - %6.4f\n"%(n_omitted_cl,n_omitted_confs, first_omitted_e, last_omitted_e) + clu_str
        return clu_str, dlg_run_str



    def read(self, filename):
        """
1.5 2.0 binding
  0   0   0   0  0.000  0.000  0.000 -14.645
  1   0   1   0  2.449  0.000  0.000 -14.636
  0   1   0   1  1.281  1.281  1.281 -14.424
  1   1   1   1  2.548  1.014  1.014 -14.210
                   |      |      |
                   |      |      |-rmsd from cluster seed @ 2.0 tolerance
                   |      |-rmsd from cluster seed @ 1.0 tolerance
                   |-rmsd from overall reference structure      
        """
        file_ptr = open(filename)
        lines = file_ptr.readlines() # read the file
        file_ptr.close()

        d = self.clustering_dict # local copy

        # see if there's the last char is a 'd' or 'b'  or 'e'to denote energy
        word_list = string.split(lines[0])
        if word_list[-1][0] in ['b', 'd', 'e']: 
            file_energy_used = word_list[-1]
            ind = string.find(lines[0], file_energy_used)
            lines[0] = lines[0][:ind] # strip energy symbol
        else:
            file_energy_used = 'binding'

        # ??? make sure we're consistent with the argsort
        #check that self has done some clustering before this
        if len(self.clustering_dict)>0:
            assert file_energy_used[0] == self.energy_used[0], 'Cluster energy mismatch'
        else:
            #here the dlg had no clustering in it
            self.energy_used = file_energy_used
            #redo argsort with file_energy_used
            if file_energy_used[0] == 'd':
                energy_list = [conf.docking_energy for conf in self.data]
            else:
                energy_list = [conf.binding_energy for conf in self.data]
            # sort the conformations by energy
            self.argsort = Numeric.argsort(energy_list)
            #print "self.argsort=", self.argsort

        t_list = map(float, string.split(lines[0]))
        #t_list for the example is [0.5, 2.0]
        num_t = len(t_list)   #number of clusterings

        for tolerance in t_list: # initialize the keys
            if d.has_key(tolerance):
                raise RuntimeError, "overwriting existing clustering"
            c = d[tolerance] = Clustering()
            c.tolerance = tolerance

        # cx is the index into self.data, the list of conformatons
        # NEW FORMAT: 
        # 0   0   0   0  0.000  0.000  0.000 -14.645
        #first line has list of tolerances
        for cx, l in enumerate(lines[1:]):
            ll = l.split()
            #eg: 2 tolerances gives [ 0,  0,  0,  0, 0.000, 0.000, 0.000,-14.645]
            #num_t *2
            c_list = map(int, ll[:num_t*2])
            data_list = map(float, ll[num_t*2:])

            for t, i in zip(t_list, xrange(len(c_list)/2)):
                cluster_index = c_list[2*i]
                cluster_rank = c_list[2*i+1]
                conf = self.data[int(self.argsort[cx])]
                if cluster_rank == 0:
                    assert len(d[t]) == cluster_index
                    d[t].append(Cluster(conf))
                else:
                    # add conformation to its cluster
                    assert len(d[t][cluster_index]) == cluster_rank
                    d[t][cluster_index].append(conf)
                # tell the conformation what cluster(s) it belongs to...
                conf.cluster_dict[t] = (cluster_index, cluster_rank)
                conf.refRMS = data_list[0]
                conf.clRMS = data_list[1]
                #FIX THIS: this ends up in setting it for the last clustering value only
                #should this be set here????... i think not
                #setattr(conf, self.energy_used, data_list[-1])


    def rebuild_clusters(self, clusterLists, tolerance):
        """
        clusterLists: 
            ordered list of cl_lists built by parsing dlg
            first cl_list is for lowest energy cluster
            second cl_list is for next lowest energy cluster, etc
        cl_list:
            ordered lists of info for cluster's members:
                [clrank, rank, run, energy, clRMS, refRMS]
            within each cl_list, conformations are also ranked by energies.
        """

        clusters = Clustering()
        
        for cl_list in clusterLists:
            #new Cluster is initialized with seed=conf of lowest energy
            #initialize Cluster with 1st conf which is cl_list[0]
            # then 2 is index of run info which is 1 based
            ind = cl_list[0][2]-1
            conf = self.data[ind]
            c = Cluster(conf)
            init = 1
            for l in cl_list:
                #use l[0], run number, to get index to conf
                conf = self.data[l[2]-1]
                #print 'using conf ', self.data[l[2]-1]
                if not init:
                    c.append(conf)
                init = 0
                conf.energy, conf.clRMS, conf.refRMS = l[3:]
            clusters.append(c)
        #should this be extend or add or something?
        self.clustering_dict[tolerance] = clusters
        #tell all the conformations about this_clustering: 'clusters'
        for (cx, cl) in map(None, xrange(len(clusters)),
                                 clusters):
            for (rank, conf) in map(None, xrange(len(cl)), cl):
                conf.cluster_dict[tolerance] = (cx, rank)



