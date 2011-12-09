#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2010
#
# adapted from earlier code from B. Norledge
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/sequence.py,v 1.1 2010/09/13 22:21:33 sanner Exp $
# 
# $Id: sequence.py,v 1.1 2010/09/13 22:21:33 sanner Exp $
#

import numpy
from MolKit import pdbWriter

from MolKit.PDBresidueNames import AAnames
oneLetterNames = AAnames
threeLetterNames = oneLetterNames.keys()


class Sequence:

    def __init__(self,sequence=None, numbers=None, name='None'):
        """ numbers is an optional list, same length as sequence,
        with the corresponding
        sequence numbers. Note that gaps are also numbered
        """
        self.name = name
        self.sequence = []
        if sequence:
            for i in range(len(sequence)):
                resName = sequence[i]
                if len(resName)!=1:
                    if resName in threeLetterNames:
                        residue = oneLetterNames[resName]
                    elif '-' in resName:
                        residue = '-'
                    else:
                        residue = 'X'
                else:
                    residue = resName.upper()
                self.sequence.append(residue)
        self.applyNumbers(numbers)


    def applyNumbers(self, numbers=None):
        gapMap = map(lambda x: x in ['-','|'],self.sequence)
        ngaps = numpy.sum(gapMap)
        nresidues = len(gapMap)-ngaps
        if numbers is None:
            numbers = map(lambda x: str(x+1),range(nresidues))
        if len(numbers)!=nresidues:
            raise ValueError('Numbers do not correspond to all residues')
        self.numbers = numbers
        count=0
        newnumbers = []
        for i in gapMap:
            if i:
                newnumbers.append('')
            else:
                newnumbers.append(str(numbers[count]))
                count=count+1
        self.gappednumbers = newnumbers
    
    def __repr__(self):
        repr = ''
        for residue in self.sequence[:10]:
            repr = repr+residue
        repr = '<Sequence instance> %s: %10s...' % (self.name,repr)
        return repr

    def __len__(self):
        return len(self.sequence)
    
    def __add__(self,other):
        """ Currently this will renumber everything from scratch to avoid duplication of residue numbers
        """
        sequence = self.sequence + other.sequence
        return Sequence(name=self.name,sequence=sequence)

    def __getitem__(self,index):
        return self.sequence[index]



class Alignment:
    """ Base class for a sequence alignment. Data is a dictionary with molecule
    identifiers as keys and the aligned sequences as values"""

    def __init__(self, sequences=None, name=None):
        self.name = name
        self.sequences={}
        self.seqNames=[]
        if sequences:
            for sequence in sequences:
                self.addSequence(sequence)
        self.writer = pdbWriter.PdbWriter()
        
    def __repr__(self):
        if len(self.sequences)>0:
            repr =  '<Alignment instance> with %d sequences of length %d:' % (
                len(self.sequences),len(self))
        else:
            repr =  '<Alignment instance> with 0 sequences'
        return repr

    def __len__(self):
        if len(self.sequences)>0:
            return len(self[0])
        return 0

    def __add__(self,other):
        new_aln = Alignment()
        for seqName in self.seqNames:
            new_aln.addSequence(self.sequences[seqName])
        for seqName in other.seqNames:
            new_aln.addSequence(other.sequences[seqName])
        return new_aln
    
    def __getitem__(self,index):
        if type(index)==type(1):
            #return a single sequence:
            seqName = self.seqNames[index]
            return self.sequences[seqName]
        elif type(index)==type(slice(1)):
            #return another alignment
            aln = Alignment()
            seqNames = self.seqNames[index.start:index.stop]
            for seqName in seqNames:
                aln.addSequence(sequence=self.sequences[seqName])
            return aln
        
    def read(self,alnFileName):
        data = open(alnFileName).readlines()
        if data[0][:7]!='CLUSTAL':
            print 'Not a clustalformatted file'
            return None
        sequences = {}
        seqNames = []
        for line in data[1:]:
            if line[0].isalnum():
                info = line.split()
                seqName = info[0]
                seqData = info[1]
                if not sequences.has_key(seqName):
                    sequences[seqName]=Sequence(name=seqName)
                    seqNames.append(seqName)
                sequences[seqName] = sequences[seqName]+Sequence(sequence=seqData)
        for seqName in seqNames:
            #if seqName matches a current sequence, need to rehash this sequence with the new
            # arrangement of gaps, but keep the old numbers
            #NB -this assumes the old and new sequences have the same number of residues!!
            index = None
            if seqName in self.seqNames:
                # get the old numbering, stripped of gaps
                sequence = self.sequences[seqName]
                numbers = []
                for number in sequence.numbers:
                    if number != '': numbers.append(number)
                #remove the old copy of the sequence (with the old gaps) from the alignment
                index = self.seqNames.index(seqName)
                self.deleteSequence(seqName)
                #replace the read sequence's numbering system with the new one
                sequences[seqName].applyNumbers(numbers)
            self.addSequence(sequences[seqName],index)

    def write(self,alnFileName):
        outfile = open(alnFileName,'w')
        title = 'CLUSTAL W multiple sequence alignment\n\n'
        outfile.write(title)
        nsegments = int(len(self)/60.)+1
        for x in range(nsegments):
            for sequence in self:
                outstring = sequence.name.ljust(16)
                for residue in sequence.sequence[60*x:60*x+60]:
                    outstring = outstring + residue
                outfile.write(outstring+'\n')
            outfile.write('\n\n')
        outfile.close()
        
    def trim(self):
        """get rid of any universal gaps in the alignment.
        """
        #make sure we have an alignment
        if len(self)==0:
            return
        #make sure we have an up-to-date matrix
        if not hasattr(self,'matrix'):
            self.makeMatrix()
        nsequences,nresidues = numpy.shape(self.matrix)
        if (nsequences != len(self.sequences) or
            nresidues != len(self)):
            self.makeMatrix()
        transpose = numpy.transpose(self.matrix)
        gaplist = []
        #any row with sum=0 in the transpose corresponds to a column in the alignment
        #which is all gaps. So add the positions of these columns to the gaplist
        for x in range(len(transpose)):
            line = transpose[x]
            if numpy.sum(line)==0:
                gaplist.append(x)
        #now can simply pop the unwanted gaps out of each sequence.
        gaplist.reverse()            
        for sequence in self:
            for gap in gaplist:
                junk=sequence.sequence.pop(gap)
                junk=sequence.gappednumbers.pop(gap)

    def makeMatrix(self):
        """ Sets up a matrix (nsequences x len(sequences), the elements of which are 0 for a gap,
        1 for anything else. Used by the trim command"""
        self.matrix = []
        for x in range(len(self.sequences)):
            numbers = self[x].gappednumbers
            line = map(lambda x: x != '',numbers)
            self.matrix.append(line)
        self.matrix = numpy.array(self.matrix)

    
    def addSequence(self, sequence, index=None):
        """add a sequence to the alignment. Gets tagged on to the end
        unless index is supplied, when it will be inserted at that
        position
        """
        if self.sequences:
            difflen = len(sequence)-len(self)
            if difflen >0:
                seqNames = self.sequences.keys()
                addOn = difflen*'-'
                for seqName in seqNames:
                    self.sequences[seqName] = self.sequences[seqName]+Sequence(sequence=addOn)
            elif difflen <0:
                difflen = -difflen
                addOn = difflen*'-'
                sequence = sequence + Sequence(sequence=addOn)
        self.sequences[sequence.name]=sequence
        if index is None:
            self.seqNames.append(sequence.name)
        else:
            self.seqNames = self.seqNames[:index]+[sequence.name]+self.seqNames[index:]

    def deleteSequence(self,sequenceName):
        if sequenceName not in self.seqNames:
            return
        del(self.sequences[sequenceName])
        idx = self.seqNames.index(sequenceName)
        junk = self.seqNames.pop(idx)


