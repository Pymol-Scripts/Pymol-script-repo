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
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/seqViewer.py,v 1.3 2010/09/13 23:33:41 sanner Exp $
# 
# $Id: seqViewer.py,v 1.3 2010/09/13 23:33:41 sanner Exp $
#

import Pmw, Tkinter, tkFileDialog
from MolKit.sequence import Sequence, Alignment

residueColors = {'A':'black', 'C':'black', 'D':'black', 'E':'black',
                 'F':'black', 'G':'black', 'H':'black', 'I':'black',
                 'K':'black', 'L':'black', 'M':'black', 'N':'black',
                 'P':'black', 'Q':'black', 'R':'black', 'S':'black',
                 'T':'black', 'V':'black', 'W':'black', 'Y':'black',
                 '?':'green'}


class AlignmentEditor(Tkinter.Frame):
    """
    GUI for editing sequence alignments. Note to self (and anyone
    else who cares...): the top thing on the window is bottom of the
    displayList..."""

    def __init__(self, alignment=None, master=None, name=None, **kw):
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
            master = Tkinter.Toplevel(self.Master)
            master.title('Alignment Editor')
        else:
        #    master = Tkinter.Toplevel(self.Master)
            master = self.Master
        Tkinter.Frame.__init__(self ,master)
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
        #print 'Filling Canvas'
        yCoord = 0
        sequences = self.alignment.sequences
        seqNames = self.alignment.seqNames
        seqCount=0
        self.gapnum=gapnum=0
        for seqName in seqNames:
            seqCount = seqCount+1
            seqTag = seqName
            sequence = sequences[seqName].sequence
            numbers = sequences[seqName].gappednumbers
            resTag = 'name'
            #print seqTag,resTag
            uniqTag = seqTag+'_'+resTag
            self.canvas.create_text(0, yCoord, text=seqName,
                                    tags=(resTag,seqTag,uniqTag))
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
