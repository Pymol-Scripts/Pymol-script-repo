## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/alignmentEditor.py,v 1.7 2007/07/24 17:30:40 vareille Exp $
#
# $Id: alignmentEditor.py,v 1.7 2007/07/24 17:30:40 vareille Exp $
#



import os
import numpy.oldnumeric as Numeric
from mglutil.math import rigidFit
from MolKit import pdbWriter
import Tkinter
import tkFileDialog
import string
import Pmw


oneLetterNames = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
                  'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
                  'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
                  'TRP':'W','TYR':'Y'}
threeLetterNames = oneLetterNames.keys()
residueColors = {'A':'black', 'C':'black', 'D':'black', 'E':'black',
                 'F':'black', 'G':'black', 'H':'black', 'I':'black',
                 'K':'black', 'L':'black', 'M':'black', 'N':'black',
                 'P':'black', 'Q':'black', 'R':'black', 'S':'black',
                 'T':'black', 'V':'black', 'W':'black', 'Y':'black'}

Fitter = rigidFit.RigidfitBodyAligner


class Sequence:

    def __init__(self,sequence=None, numbers=None, name=None):
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
                    residue = string.upper(resName)
                self.sequence.append(residue)
        self.applyNumbers(numbers)

    def applyNumbers(self,numbers=None):
        gapMap = map(lambda x: x in ['-','|'],self.sequence)
        ngaps = Numeric.sum(gapMap)
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

    def __init__(self,sequences=None,name=None):
        self.name = name
        self.sequences={}
        self.seqNames=[]
        if sequences:
            for sequence in sequences:
                self.addSequence(name=sequence.name,sequence=sequence)
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
        nsequences,nresidues = Numeric.shape(self.matrix)
        if (nsequences != len(self.sequences) or
            nresidues != len(self)):
            self.makeMatrix()
        transpose = Numeric.transpose(self.matrix)
        gaplist = []
        #any row with sum=0 in the transpose corresponds to a column in the alignment
        #which is all gaps. So add the positions of these columns to the gaplist
        for x in range(len(transpose)):
            line = transpose[x]
            if Numeric.sum(line)==0:
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
        self.matrix = Numeric.array(self.matrix)

    
    def addSequence(self,sequence,index=None):
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

class AlignmentEditor(Tkinter.Frame):
    """ GUI for editing sequence alignments. Note to self (and anyone
    else who cares...): the top thing on the window is bottom of the
    displayList..."""

    def __init__(self,alignment=None,master=None,name=None,**kw):
        self.name=name
        self.xspace = 10
        self.yspace = 20
        self.selection = ()
        self.selectionTags = []
        self.colors = residueColors
        self.colors['default']='black',
        self.colors['selection']='yellow'
        self.colors['|']='magenta'
        #somewhere to store any tagspecifc colors
        self.colors['special']={}
        self.Master=master
        if alignment:
            self.alignment = alignment
            self.name = alignment.name
        else:
            self.alignment = Alignment()
        self.createGUI()

    def createGUI(self):
        #self.widgetArea = Tkinter.Frame(self, borderwidth=2, relief='sunken')
        #if self.hasGUI:
        #    return
        if self.Master is None:
            master = Tkinter.Tk()
        else:
            master = Tkinter.Toplevel(self.Master)
        master.title('Alignment Editor')
        Tkinter.Frame.__init__(self,master)
        Tkinter.Pack.config(self, expand=1, fill=Tkinter.BOTH)
        self.createMenus()
        self.canvasFrame = Tkinter.Frame(self)
        self.canvas = Pmw.ScrolledCanvas(self.canvasFrame,usehullsize=1,
                                         hull_width=600,hull_height=200,
                                         hscrollmode='dynamic',
                                         vscrollmode='dynamic',
                                         canvasmargin=1)
        self.canvas.pack(side=Tkinter.LEFT,expand=1,fill=Tkinter.BOTH)
        self.canvasFrame.pack(side=Tkinter.LEFT,expand=1,fill=Tkinter.BOTH)
        self.canvas._canvas.bind("<ButtonPress-1>",self.mouseDown)
        self.canvas._canvas.bind("<Button1-Motion>",self.mouseMotion)
        self.canvas._canvas.bind("<Button1-ButtonRelease>",self.mouseUp)
        self.canvas._canvas.bind("<Shift-ButtonPress-1>",self.startSelection)
        self.canvas._canvas.bind("<Shift-Button1-Motion>",self.continueSelection)
        self.canvas._canvas.bind("<Shift-Button1-ButtonRelease>",self.mouseSelect)
        self.canvas._canvas.bind("<Control-ButtonPress-1>",self.startSelection)
        self.canvas._canvas.bind("<Control-Button1-Motion>",self.continueSelection)
        self.canvas._canvas.bind("<Control-Button1-ButtonRelease>",self.mouseDeselect)
        self.fillCanvas()
        #self.hasGUI = 1
        
    def startSelection(self,event=None):
        #print 'In startSelection'
        self.x0 = self.canvas.canvasx(event.x)
        self.y0 = self.canvas.canvasy(event.y)
        
    def continueSelection(self,event=None):
        #print 'In continueSelection'
        self.x1=self.canvas.canvasx(event.x)
        self.y1=self.canvas.canvasy(event.y)
        self.clearSelBox()
        self.canvas.create_line(self.x0,self.y0,self.x1,self.y0,tags=('selBox'))
        self.canvas.create_line(self.x0,self.y0,self.x0,self.y1,tags=('selBox'))
        self.canvas.create_line(self.x0,self.y1,self.x1,self.y1,tags=('selBox'))        
        self.canvas.create_line(self.x1,self.y0,self.x1,self.y1,tags=('selBox'))

    def clearSelBox(self):
        #print 'In clearSelBox'
        items = self.canvas.find_withtag('selBox')
        for item in items:
            self.canvas.delete(item)

    def mouseSelect(self,event=None,deselect=0):
        """ deselect=1 if removing selection
        """
        #print 'In select'
        self.x1=self.canvas.canvasx(event.x)
        self.y1=self.canvas.canvasy(event.y)
        self.clearSelBox()
        items = self.canvas.find_overlapping(self.x0,self.y0,self.x1,self.y1)
        #print items
        self.select(items,deselect)

    def select(self,items,deselect=0):
        if not items:
            return
        self.lastSelect = [deselect,[]]
        for item in items:
            resTag,seqTag,uniqTag = self.canvas.gettags(item)[:3]
            if not deselect:
                if item not in self.selection:
                    self.canvas.itemconfig(item,fill=self.colors['selection'])
                    self.canvas.addtag_withtag('selected',uniqTag)
            else:
                if item in self.selection:
                    self.canvas.itemconfig(item,fill=self.colors['default'])
                    self.canvas.dtag(item,'selected')
            self.lastSelect[1].append(item)
        self.selection = self.canvas.find_withtag('selected')
        self.rebuildSelectionTags()

        
    def rebuildSelectionTags(self):
        #print 'In rebuildSelectionTags'
        self.selectionTags = []
        for item in self.selection:
            self.selectionTags.append(self.canvas.gettags(item))

            
    def mouseDeselect(self,event=None):
        #print 'In deselect'
        self.mouseSelect(event=event,deselect=1)
        

    def clearSelection(self):
        #print 'In clearSelection'
        self.updateColor(self.selectionTags,self.colors['default'])
        self.canvas.dtag('all','selected')
        self.selection = ()
        self.selectionTags = []
        

    def mouseDown(self, event=None):
        #print 'In mouseDown'
        # this method has to figure out where we are on the canvas when the button is pressed
        tags = self.canvas.gettags('current')
        #markers for the mousemotion
        self.x0 = self.canvas.canvasx(event.x)
        self.y0 = self.canvas.canvasy(event.y)
        #return
        if tags:
            self.resTag = tags[0]
            self.seqTag = tags[1]
            self.uniqTag = tags[2]
            self.currentResidue=self.canvas.find_withtag(self.uniqTag)
            self.currentSeq=self.canvas.find_withtag(self.seqTag)
            self.findAllToRight()
            self.findToLeft()
            self.findNeighborSequences()


    def findNeighborSequences(self):
        currTag = self.seqTag
        currIndex = self.alignment.seqNames.index(currTag)
        if currIndex != len(self.alignment.seqNames)-1:
            nextTag = self.alignment.seqNames[currIndex+1]
            self.nextSeq = self.canvas.find_withtag(nextTag)
        else:
            self.nextSeq = ()
        if currIndex != 0:
            prevTag = self.alignment.seqNames[currIndex-1]
            self.prevSeq = self.canvas.find_withtag(prevTag)
        else:
            self.prevSeq = ()

        
    def findAllToRight(self):
        #print 'In findAllToRight'
        currentSeq = list(self.currentSeq)
        index= currentSeq.index(self.currentResidue[0])
        self.allToRight = tuple(currentSeq[index:])

        
    def findToLeft(self):
        #print 'In findToLeft'
        currentResidue = self.currentResidue
        if currentResidue[0] == self.canvas.find_all()[0]:
            return None
        prevResidue = self.canvas.find_below(currentResidue)
        tags = self.canvas.gettags(prevResidue)
        self.toLeft = prevResidue


    def mouseMotion(self, event=None):
        #print 'In mouseMotion'
        self.x1 = self.canvas.canvasx(event.x)
        self.y1 = self.canvas.canvasy(event.y)
        tags = self.canvas.gettags(self.currentResidue)
        #dragging to left closes a gap. Can't be done if there is no gap.
        if (self.x1 < self.x0-self.xspace and
            self.canvas.itemcget(self.toLeft,'text')=='-' and
            'movable' in tags):
            self.closeGap()
        #dragging to right opens a gap. Can always be done.
        if (self.x1 > self.x0+self.xspace and
            'movable' in tags):
            self.openGap()
        #dragging up swaps current and previous sequences
        if self.y1 > self.y0+self.yspace:
            self.swapSequences(self.currentSeq,self.nextSeq,direction='down')
        if self.y1 < self.y0-self.yspace:
            self.swapSequences(self.prevSeq,self.currentSeq,direction='up')
        
        
    def swapSequences(self,topSequence,bottomSequence,direction):
        #print 'In swapSequences'
        #check both sequences exist:
        if not (topSequence and bottomSequence):
            return
        currTag = self.seqTag
        seqIndex = self.alignment.seqNames.index(currTag)
        #move the sequences up and down:
        for item in topSequence:
            self.canvas.move(item,0,self.yspace)
        for item in bottomSequence:
            self.canvas.move(item,0,-self.yspace)
        #reset the starting coordinates and update the
        #displaylist. Also update the order of the alignment sequences
        if direction=='down':
            #swapping currentSeq and nextSeq
            nextTag = self.alignment.seqNames[seqIndex+1]
            self.alignment.seqNames[seqIndex]=nextTag
            self.alignment.seqNames[seqIndex+1]=currTag
            self.y0 = self.y0+self.yspace
            items = list(self.currentSeq)
            items.reverse()
            for item in items:
                self.canvas._canvas.lift(item,(self.nextSeq[-1],))
        else:
            #swapping prevSeq and currentSeq
            prevTag = self.alignment.seqNames[seqIndex-1]
            self.alignment.seqNames[seqIndex]=prevTag
            self.alignment.seqNames[seqIndex-1]=currTag
            self.y0 = self.y0-self.yspace
            items = list(self.prevSeq)
            items.reverse()
            for item in items:
                self.canvas._canvas.lift(item,(self.currentSeq[-1],))
        self.findToLeft()
        self.findNeighborSequences()

        
    def closeGap(self):
        # delete the gap item
        self.canvas.delete(self.toLeft)
        # update self.toLeft
        self.findToLeft()
        # move current item to left
        for item in self.allToRight:
            self.canvas.move(item, -self.xspace, 0)
        #need to tag this sequence as edited
        tags = self.canvas.gettags(self.currentResidue)
        if 'edited' not in tags:
            self.canvas.addtag_withtag('edited',tags[1])        
        #then need to reset x0, so it can be repeated:
        self.x0 = self.x0-self.xspace
        
    def openGap(self):
        #find out where we are so we know where to insert
        coordx,coordy = self.canvas.coords(self.currentResidue)
        #move everything that is to the right, even further to the right
        for item in self.allToRight:
            self.canvas.move(item, self.xspace, 0)
        #insert a gap at the current coordinates
        tags = self.canvas.gettags(self.currentResidue)
        restag,seqtag,uniqtag = tags[:3]
        restag = 'gap'+str(self.gapnum)
        self.gapnum = self.gapnum+1
        uniqtag = seqtag+'_'+restag
        #need to tag this particular sequence as edited
        if 'edited' not in tags:
            self.canvas.addtag_withtag('edited',seqtag)
        #need to add the new gap, also tagged as edited
        self.canvas.create_text(coordx,coordy,
                                text='-',
                                tags=(restag,seqtag,uniqtag,'movable','edited'))
        #need to resize scroll region to accomodate new gap
        self.canvas.resizescrollregion()
        #then need to move things around in the displaylist so that the new item
        #is just below the old one (i.e just above toLeft)
        newItem = self.canvas.find_all()[-1]
        self.canvas._canvas.lift(newItem,self.toLeft)
        #then need to redefine the current sequence
        self.currentSeq=self.canvas.find_withtag(self.seqTag)
        #need to find the new toLeft
        self.findToLeft()
        #then need to reset x0, so it can be repeated:
        self.x0 = self.x0+self.xspace   

    def mouseUp(self, event=None):
        self.remakeAlignment()
        #self.redraw()
        return


    def remakeAlignment(self):
        """Replaces the edited sequences in the underlying alignment.
        """
        seqStr = None
        edited = self.canvas.find_withtag('edited')
        if edited ==():
            return
        for item in edited:
            tags = self.canvas.gettags(item)
            if tags[0] == 'name':
                if seqStr:
                    sequence = Sequence(name=sequenceName,sequence=seqStr)
                    self.alignment.deleteSequence(sequenceName)
                    self.alignment.addSequence(sequence,index)
                sequenceName = tags[1]
                index = self.alignment.seqNames.index(sequenceName)
                seqStr=''
            else:
                seqStr = seqStr+self.canvas.itemcget(item,'text')
        sequence = Sequence(name=sequenceName,sequence=seqStr) #tag on the final sequence
        self.alignment.deleteSequence(sequenceName)
        self.alignment.addSequence(sequence,index)
        self.canvas.dtag('all','edited')
        
    def createMenus(self):
        #print 'In createMenus'
        self.mBar = Tkinter.Frame(self, relief=Tkinter.RAISED,borderwidth=2)
        self.mBar.pack(fill=Tkinter.X)
        self.menuButtons = {}
        self.makeFileMenu()
        self.makeEditMenu()
        apply(self.mBar.tk_menuBar, self.menuButtons.values())
        self.title = Tkinter.Label(self.mBar, text=self.name)
        self.title.pack(side=Tkinter.RIGHT)
        
    def makeFileMenu(self):
        #print 'In makeFileMenu'
        File_button = Tkinter.Menubutton(self.mBar, text='File',underline=0)
        self.menuButtons['File'] = File_button
        File_button.pack(side = Tkinter.LEFT, padx='1m')
        File_button.menu = Tkinter.Menu(File_button)
        File_button.menu.add_command(label='Load...', underline=0,
                                     command = self.loadFile)
        File_button.menu.add_command(label='Write...', underline=0,
                                     command = self.writeFile)
        File_button.menu.add_command(label='Exit...', underline=0,
                                     command = self.exit)
        File_button['menu'] = File_button.menu

    def loadFile(self):
        #print 'In loadFile'
        title = 'Read CLUSTAL formatted alignment file'
        types = [('CLUSTAL files', '*.aln')]
        file = tkFileDialog.askopenfilename( filetypes=types,
                                             title=title)
        if file:
            self.alignment.read(file)
        self.redraw()


    def writeFile(self):
        #print 'In writeFile'
        #self.remakeAlignment() # always done in mouseUp
        #self.redraw() # horribly expensive
        title = 'Save CLUSTAL formatted alignment file'
        types = [('CLUSTAL files', '*.aln')]
        file = tkFileDialog.asksaveasfilename( filetypes=types,
                                               title=title)
        if file and self.alignment:
            self.alignment.write(file)


    def fillCanvas(self):
        print 'Filling Canvas'
        yCoord = 0
        sequences = self.alignment.sequences
        seqNames = self.alignment.seqNames
        seqCount=0
        self.gapnum=gapnum=0
        for seqName in seqNames:
            seqCount=seqCount+1
            seqTag = seqName
            sequence=sequences[seqName].sequence
            numbers =sequences[seqName].gappednumbers
            resTag = 'name'
            #print seqTag,resTag
            uniqTag = seqTag+'_'+resTag
            self.canvas.create_text(0,yCoord,text=seqName,tags=(resTag,seqTag,uniqTag))
            for xCoord in range(len(self.alignment)):
                resName = sequence[xCoord]
                resTag = numbers[xCoord]
                #need unique tags for gaps too
                if resTag == '':
                    resTag = 'gap'+str(gapnum)
                    gapnum=gapnum+1
                uniqtag = seqTag+'_'+resTag
                try:
                    fillColor = self.colors[resName]
                except:
                    fillColor = self.colors['default']
                self.canvas.create_text(100+xCoord*self.xspace,yCoord,
                                        text=resName,
                                        fill=fillColor,
                                        tags=(resTag,seqTag,uniqtag,'movable'))
            yCoord = yCoord+self.yspace
        self.canvas.resizescrollregion()
        print 'updating colors'
        self.updateColor(self.selectionTags,'yellow')
        self.updateSpecialColor()
        print 'Done'
        
    def updateSpecialColor(self):
        if self.colors['special']=={}:
            return
        for tag in self.colors['special'].keys():
            item = self.canvas.find_withtag(tag)[0]
            self.canvas.itemconfig(item,fill=self.colors['special'][tag])

    def updateColor(self,tags,color):
        #print 'In updateColor'
        if not tags:
            return
        #print tags
        for tag in tags:
            item = self.canvas.find_withtag(tag[2])
            self.canvas.itemconfig(item,fill=color)

    def exit(self, event=None):
        if self.Master is not None:
            self.master.withdraw()
        else:
            self.master.destroy()   
        
    def makeEditMenu(self, event=None):
        #print 'In makeEditMenu'
        Edit_button = Tkinter.Menubutton(self.mBar, text='Edit',underline=0)
        self.menuButtons['Edit'] = Edit_button
        Edit_button.pack(side = Tkinter.LEFT, padx='1m')
        Edit_button.menu = Tkinter.Menu(Edit_button)
        Edit_button.menu.add_command(label='Redraw', underline=0,
                                    command = self.redraw)
        Edit_button.menu.add_command(label='Clear Selection', underline=0,
                                    command = self.clearSelection)
        Edit_button.menu.add_command(label='Delete Selected Sequences', underline=0,
                                    command = self.deleteSelectedSequences)
        #Edit_button.menu.add_command(label='Delete Selected Residues', underline=0,
        #                            command = self.deleteSelectedResidues)
        Edit_button.menu.add_command(label='Trim Gaps', underline=0,
                                    command = self.trim)
        Edit_button['menu'] = Edit_button.menu


    def trim(self):
        self.alignment.trim()
        self.redraw()

    def deleteSelectedSequences(self):
        if self.selectionTags==[]:
            return
        while self.selection:
            firstSeqName = self.canvas.gettags(self.selection[0])[1]
            self.deleteSequence(firstSeqName)
        self.redraw()

    def deleteSequence(self,seqName):
        self.alignment.deleteSequence(seqName)
        items = self.canvas.find_withtag(seqName)
        if len(items):
            for item in items:
                tags = self.canvas.gettags(item)
                self.canvas.delete(item)
                try:
                    del self.colors['special'][tags[2]]
                except:
                    continue
            self.selection = self.canvas.find_withtag('selected')
            self.rebuildSelectionTags()
            
    def deleteSelectedResidues(self):
        #print 'In deleteSelection'
        #if selection is empty get out
        if self.selectionTags==[]:
            return
        #get hold of a list of selected residues
        selSeq = map(lambda x: x[1], self.selectionTags)
        #uniquify it
        uniqSelSeq = [selSeq[0]]
        for seq in selSeq[1:]:
            if seq != uniqSelSeq[-1]:
                uniqSelSeq.append(seq)
        #for each sequence, build a new sequence minus the selected tags, and
        #update the alignment
        for seqName in uniqSelSeq:
            sequence = []
            residues = self.canvas.find_withtag(seqName)
            for residue in residues[1:]:
                tags = self.canvas.gettags(residue)
                if 'selected' not in tags:
                    sequence.append(self.canvas.itemcget(residue,'text'))
            sequence = Sequence(name=seqName,sequence=sequence)
            index= self.alignment.seqNames.index(seqName)
            self.alignment.deleteSequence(seqName)
            self.alignment.addSequence(sequence,index)
        self.selectionTags=[] # can't have any tags if the selection is all gone...
        self.redraw()
            
    def redraw(self, event=None):
        self.canvas.delete('all')
        self.fillCanvas()

if __name__ == '__main__':
    aln = AE.Alignment()
    aln.read('pdb1tab.aln')
    edt = AE.AlignmentEditor(alignment=aln)

