## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
#$Id
#
########################################################################
#Authors: Sowjanya Karnati,Michel F Sanner
#
########################################################################
#
#
#




import numpy.oldnumeric as Numeric, math, types
from molecule import Atom, AtomSet
from protein import Residue
from types import ListType, TupleType, StringType, IntType, FloatType, LongType
from numpy.oldnumeric import ArrayType, sum
"""
    GapFinder contain class to find Gap between residues.
    findGap: This class finds the residues which are not connected.
    Finds CA atoms computes distances between CA atoms ,computes mean and
    standard deviation,max distance two CA atoms can be seperated .Compares
    all distances with that and if any distance exceeds then there is a gap.
    That chain posses attribute hasGap and atoms have gap attribute.
"""

class FindGap:

    def __init__(self,mat,chain,atoms=[]):
        """
        None<--FindGap(matrix,chain)
        matrix :transformed coords
        chain:mol.chain 
        atoms:a set of ca atoms
        This function computes distances between atoms ,sd ,max 
        distance .If distance between any two atoms are greater than max then there is a gap between that residues
        """
        self.mat=mat
        self.chain=chain
        self.mol=chain.parent
        mol=self.mol
        self.atoms=atoms
        if len(mat)>1:
            self.Gapresidues = self.Gapfunction(mat,mol,atoms)   
        elif len(mat)==1:
            print "only one atom selected"
            return

        
    def Gapfunction(self,mat,mol,atoms=[]):     
            result=[] 
            #for chain in range(0,len(mol.chains)):
            residues=[]
            if atoms!=[]:
                for at in atoms:
                    if at.name!="CA":
                        print "All atoms given are not CA atoms"
                        return
            if atoms==[]:
                atoms=self.chain.getAtoms().get(lambda x: x.name=='CA')
                #print atoms
                if len(atoms)==0 and len(mol.chains)==0:
                    self.warningMsg("No CA atoms in %s"%mol.name) 
                    return
                elif len(atoms)==0 and len(mol.chains)>1:
            
                    print "No CA atoms in %s .%s"%(mol.name,self.chain.name)
                    return    
            
            atoms.sort()
            if hasattr(self.chain,"hasGap"):
                delattr(self.chain,"hasGap")
            for at in atoms:
                if hasattr(at,"gap"):
                    delattr(at,"gap")
              
            distances = self.measure_distance(mat)
            sd,mean = self.standarddeviation_of_distances(distances)
            ats = self.compute_max_min(mean,distances,mat,mol,atoms)
            if ats==[]:
                return
            for at in ats:
                x =at[1].findType(Residue)
                residues.append(x.name)
                at[1].gap='end'
                res1=at[1].parent.name
                x =at[0].findType(Residue)
                residues.append(x.name)
                at[0].gap='start'
                res2=at[0].parent.name
                self.chain.hasGap=True
                #if len(mol.chains)>1:
                result.append([self.chain.name,res1,res2])
                #if len(mol.chains)==1:
                    #result.append([res1,res2])
            return result

            
    def measure_distance(self,mat):   
        """This function calls measure distance command and return distances
        between atoms
        """
        if len(mat)==1:
            print "chain has only one CAatom"
            return
        self.dists =[]
        for num in range(0,len(mat)):
            if num+1<=len(mat)-1:
                c1 = mat[num]
                c2 = mat[num+1]
                d = c2 - c1
                self.dists.append(math.sqrt(sum(d*d)))
        return self.dists     


    def standarddeviation_of_distances(self,distances,mean=None):   
        """This function computes standard deviations of distances between
        atoms list given
        """
        if len(distances)==1:
            mean=3.50
            self.stddev=0.2
        else:
            sum =0
            for dis in distances:
                sum =sum+dis
            #finding mean    
            mean =sum/len(distances)
            if mean>4.00:
                mean=3.50
            self.distsquares =0
            self.midval=0
            self.midval1=0
            #Now, subtract the mean individually from each of the numbers and square it 
            for dist in distances:
                self.distsquares =0
                self.distsquares=math.pow((dist-mean),2)
                self.midval=self.midval+self.distsquares
            if len(distances)>1:
                self.midval1 =    self.midval/(len(distances)-1)
            if len(distances)==1:
                mean= distances[0]
                self.stddev=0
            self.stddev=math.sqrt(self.midval1)
            
        return self.stddev,mean
       

    def compute_max_min(self,mean,distances,mat,mol,atoms=[]):
        """This function computes max and min value and finds the atoms where
        distance between them is more than max and   residues which are not
        connected
        """
        min = mean-0.300
        max = mean+0.500
        ats=[]
        if len(distances)==1:
            if distances[0]>max:
                caats=self.chain.getAtoms().get(lambda x :x.name=="CA")
                ats.append([caats[0],caats[1]])
        else:
            for d in distances:
                if d>max:
                    ind =distances.index(d)
                    if atoms==[]:
                        caats=self.chain.getAtoms().get(lambda x :x.name=="CA")
                    else:
                        caats=atoms
                    if ind+1<=len(caats)-1:
                        ats.append([caats[ind],caats[ind+1]])
        
        return ats
            
        
   
    
        





