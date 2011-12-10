# Authors: Sargis Dallakyan and Rodney Harris 
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/AutoLigandCommand.py,v 1.24.4.2 2011/01/06 21:18:25 rhuey Exp $
#
# $Id: AutoLigandCommand.py,v 1.24.4.2 2011/01/06 21:18:25 rhuey Exp $
from ViewerFramework.VFCommand import CommandGUI
from Pmv.mvCommand import MVCommand, MVAtomICOM
import Pmw, os, glob, sys, tkMessageBox, Tkinter, subprocess, tkFileDialog
from MolKit.molecule import Atom
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ExtendedSliderWidget
from autostartCommands import menuText
import AutoDockTools 
from mglutil.popen2Threads import SysCmdInThread
import cPickle
from Pmv.pmvPalettes import AtomElements
from MolKit.radii_patterns import AAradii
import copy
AutoLigandPath = os.path.join(AutoDockTools.__path__[0],'AutoLigand.py') 
from MolKit.protein import Protein, Chain, Residue
from MolKit.molecule import Atom, AtomSet
if not os.path.exists(AutoLigandPath):
    sys.stderr.write(AutoLigandPath + "not found\n")

from mglutil.gui.BasicWidgets.Tk.player import Player
from Pmv.moleculeViewer import DeleteAtomsEvent, AddAtomsEvent
from MolKit.molecule import Atom, AtomSet, Molecule
from MolKit.protein import Protein, Chain, Residue

class FloodPlayer(Player):
    def __init__(self, command, file):
        master = command.vf.GUI.ROOT
        self.autoLigandCommand = command.vf.AutoLigandCommand
        self.autoLigandCommand.spheres.Set(visible=1)
        self.autoLigandCommand.halo.Set(visible=1)
        pkl_file = open(file, 'rb')
        self.floods = []
        try:
            data = cPickle.load(pkl_file)
        except Exception, inst:
            print "Error loading ", __file__, "\n", inst
        self.xcent = data[0]
        self.ycent = data[1]
        self.zcent = data[2]
        self.centerx = data[3]
        self.centery = data[4]
        self.centerz = data[5]
        self.spacing = data[6]
        self.centers = []
        data = cPickle.load(pkl_file)
        self.floods.append(data[1])
        try:
            while data:
                data = cPickle.load(pkl_file)
                flood = copy.copy(self.floods[-1])
                for item in data[0]:
                    flood.remove(item)
                for item in data[1]:
                    flood.append(item)
                self.floods.append(flood)
        except EOFError:
            pass
        pkl_file.close()
        fileName = os.path.splitext(os.path.split(file)[-1])[0]
        self.mol = Protein(fileName)
        self.mol.allAtoms = AtomSet([])
        chain = Chain()
        self.residue = Residue(type="UNK")
        chain.adopt(self.residue, setChildrenTop=1)
        self.mol.adopt(chain, setChildrenTop=1)
        self.mol.parser = None
        self.filename = file
        fl = self.floods[0][0]
        x = (fl[1] - self.xcent)*self.spacing + self.centerx
        y = (fl[2] - self.ycent)*self.spacing + self.centery
        z = (fl[3] - self.zcent)*self.spacing + self.centerz
        if fl[4] == 7:
            atomchr = 'P'
            # note, this will color the NA atom pink (the PDB color for Phosphorus)
            radius = AAradii[13][0]
        if fl[4] == 6:
            atomchr = 'S'
            radius = AAradii[13][0]
        if fl[4] == 5:
            atomchr = 'A'
            radius = AAradii[10][0]
        if fl[4] == 4:
            atomchr = 'O'
            radius = AAradii[1][0]
        if fl[4] == 3:
            atomchr = 'N'
            radius = AAradii[4][0]
        if fl[4] == 2:
            atomchr = 'C'
            radius = AAradii[10][0]
        if fl[4] == 1:
            atomchr = 'H'      
            radius = AAradii[15][0]
        a = Atom(atomchr, self.residue, atomchr, top=self.mol)
        a._coords = [[x,y,z]]
        a._charges = {}
        a.hetatm = 1
        a.number = 0
        a.radius = radius
        self.mol.allAtoms = self.residue.atoms
        self.mol = self.autoLigandCommand.vf.addMolecule(self.mol, False)
        self.mol.levels = [Protein, Chain, Residue, Atom]
        self.autoLigandCommand.vf.displayCPK(self.mol, scaleFactor=0.4)
        self.autoLigandCommand.vf.colorByAtomType(self.mol, ['cpk'], log=0)
        self.autoLigandCommand.vf.displayLines(self.mol, negate=True, displayBO=False, lineWidth=2, log=0, only=False)
        self.colorKeys = a.colors.keys()
        maxLen =len(self.floods)-1
        Player.__init__(self, master=master, endFrame=maxLen, maxFrame=maxLen, 
                        titleStr="AutoLigand Flood Player", hasSlider=True)
        try:# withdrew SetAnim button
            self.form.ifd.entryByName['setanimB']['widget'].grid_forget()
            self.form.autoSize()
        except:
            pass
        self.nextFrame(0)
        self.form.root.protocol('WM_DELETE_WINDOW', self.hide_cb)

    def nextFrame(self, id):
        #Player.nextFrame(self, id)
        id = int(id)
        if id == self.currentFrameIndex: return
        if self.hasCounter and self.gui:
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, str(id))
            if self.hasSlider:
                self.form.ifd.entryByName['slider']['widget'].set(id)
        self.currentFrameIndex = int(id)        
        removeAtoms = AtomSet([])
        addAtoms = AtomSet([])
        
        id = int(id)
        flood = self.floods[id]
        centers = []
        materials = []
        radii = []
        prev_coords = self.mol.allAtoms.coords
        lenAtoms = len(prev_coords)
        #self.residue.atoms = AtomSet([])
        index = 0        
        #h = self.hp.heap()
        #print h
        for fl in flood:  
            x = (fl[1] - self.xcent)*self.spacing + self.centerx
            y = (fl[2] - self.ycent)*self.spacing + self.centery
            z = (fl[3] - self.zcent)*self.spacing + self.centerz
            if fl[4] == 7:
                atomchr = 'P'
                # note, this will color the NA atom pink (the PDB color for Phosphorus)
                radius = AAradii[13][0]
            if fl[4] == 6:
                atomchr = 'S'
                radius = AAradii[13][0]
            if fl[4] == 5:
                atomchr = 'A'
                radius = AAradii[10][0]
            if fl[4] == 4:
                atomchr = 'O'
                radius = AAradii[1][0]
            if fl[4] == 3:
                atomchr = 'N'
                radius = AAradii[4][0]
            if fl[4] == 2:
                atomchr = 'C'
                radius = AAradii[10][0]
            if fl[4] == 1:
                atomchr = 'H'      
                radius = AAradii[15][0]
            if not [x,y,z] in prev_coords:
                a = Atom(atomchr, self.residue, atomchr, top=self.mol)
                a._coords = [[x,y,z]]
                a._charges = {}
                a.hetatm = 1
                a.radius = radius
                #a.number = lenAtoms + 1
                addAtoms.append(a)
                lenAtoms += 1
                for key in self.colorKeys:
                    a.colors[key]=AtomElements[atomchr]
                    a.opacities[key]=1.0
            else:
                centers.append([x,y,z])
                            
#            a = Atom(atomchr, self.residue, atomchr, top=self.mol)
#            a._coords = [[x,y,z]]
#            a._charges = {}
#            a.hetatm = 1
#            a.number = index 
#            index += 1
            #aterials.append(AtomElements[atomchr])
            #enters.append([x,y,z])
            #adii.append(radius)
        #self.mol.allAtoms = self.residue.atoms
        #self.mol.geomContainer.geoms['lines'].protected = False
        #for com in self.autoLigandCommand.vf.cmdsWithOnAddObj:
        #    com.onAddObjectToViewer(self.mol)        
        #self.autoLigandCommand.vf.displayCPK(self.mol, scaleFactor=0.1)
        
        halo_centers = []
        for coord in prev_coords:
            if not coord in centers:
                index = prev_coords.index(coord)
                removeAtoms.append(self.mol.allAtoms[index])
        
        
        self.residue.assignUniqIndex() #this is needed to avoid Traceback later on
        self.mol.allAtoms.stringRepr = None #stringRepr can be very large aousing memory errors
        event = AddAtomsEvent(objects=addAtoms)
        #self.autoLigandCommand.vf.dispatchEvent(event)
        self.autoLigandCommand.vf.displayCPK.updateGeom(event)        
        event = DeleteAtomsEvent(objects=removeAtoms)
        #self.autoLigandCommand.vf.dispatchEvent(event)
        self.autoLigandCommand.vf.displayCPK.updateGeom(event)
        for atom in removeAtoms:
            self.residue.atoms.remove(atom)
        if id == self.maxFrame:
            self.autoLigandCommand.halo.Set(visible=0)
        else:
            self.autoLigandCommand.halo.Set(centers=addAtoms.coords, materials=((1,1,0,0.5),), radii=0.4)

        #self.mol.allAtoms = self.residue.atoms

        #self.vf.GUI.VIEWER.Redraw()
        #self.vf.GUI.ROOT.update()        

    def hide_cb(self):
        self.autoLigandCommand.hideGeoms()
        self.form.destroy()

        
class AutoLigandCommand(MVCommand, MVAtomICOM):        
    "GUI for AutoLigand: extends MVCommand, overwrites guiCallback"       
    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        MVAtomICOM.__init__(self)
        self.save = None
        self.showPlayer = False                    
        self.floodFile = None
        
    def onAddCmdToViewer(self):
        from DejaVu.Points import CrossSet
        self.startCross = CrossSet('StartCross', materials=((1.,1.,0),),
                            inheritMaterial=0, protected=True,
                    offset=1.0,lineWidth=5, visible=0, pickable=0,  listed=0)
        self.endCross = CrossSet('EndCross', materials=((0,1,1),),
                            inheritMaterial=0, protected=True,
                    offset=1.0,lineWidth=5, visible=0, pickable=0,  listed=0)

        from DejaVu.IndexedPolygons import IndexedPolygons
        from DejaVu.Box import Box
        from DejaVu.Spheres import Spheres
        from DejaVu import viewerConst
        from DejaVu.bitPatterns import patternList
        from opengltk.OpenGL import GL
        
        face=((0,3,2,1),(3,7,6,2),(7,4,5,6),(0,1,5,4),(1,2,6,5),(0,4,7,3))
        coords=((1,1,-1),(-1,1,-1),(-1,-1,-1),(1,-1,-1),(1,1,1),(-1,1,1),(-1,-1,1),(1,-1,1))
        #new style RGB->
        materials=((0,0,1),(0,1,0),(0,0,1),(0,1,0),(1,0,0),(1,0,0),)
        box=IndexedPolygons('Box', materials=materials, vertices=coords, faces=face,
                            inheritMaterial=0, visible=0, protected=True,  listed=0)
        box.Set(frontPolyMode=GL.GL_LINE)
        box.Set(backPolyMode=GL.GL_LINE)
        box.Set(culling=GL.GL_NONE)
        box.inheritShading=0
        box.shading=GL.GL_FLAT
        box.Set(matBind=viewerConst.PER_PART)
        box.polygonstipple.Set(pattern=patternList[0])
        box.Set(stipplePolygons=1)
        box.transparent=0
        self.box = box
        self.spheres = Spheres('Spheres', visible=0, inheritMaterial=0, radii=(0.3,), protected=True,  listed=0)
        self.halo = Spheres('Halo', visible=0, inheritMaterial=0, radii=(0.5,), protected=True,  listed=0)

                
        from DejaVu.Geom import Geom
        AutoLigand_geoms = Geom("AutoLigand_geoms", shape=(0,0),  listed=0)
        self.vf.GUI.VIEWER.AddObject(AutoLigand_geoms, parent =self.vf.GUI.miscGeom)
        self.vf.GUI.VIEWER.AddObject(self.startCross, parent=AutoLigand_geoms)
        self.vf.GUI.VIEWER.AddObject(self.endCross, parent=AutoLigand_geoms)        
        self.vf.GUI.VIEWER.AddObject(self.box, parent=AutoLigand_geoms)    
        self.vf.GUI.VIEWER.AddObject(self.spheres, parent=AutoLigand_geoms)    
        self.vf.GUI.VIEWER.AddObject(self.halo, parent=AutoLigand_geoms)    
        self.grids = {}
        self.vf.showCitation.citations["AutoLigand"] = """
Please acknowledge the use of the AutoLigand that results
in any published work, including scientific papers, films
and videotapes, by citing the following reference:
Harris R, Olson AJ, Goodsell DS.  Proteins. (2008) 70 1506..
"""
        
    def guiCallback(self, event=None):
        fileList = [] 
        fld_list = glob.glob('*.maps.fld')
        if not fld_list:
            tkMessageBox.showinfo("AutoLigand Info", "AutoLigand requires input AutoGrid maps. \nPlease click OK to select directory containing grid maps.")
            folder = tkFileDialog.askdirectory(title="Select A Folder")
            if folder:
                os.chdir(folder)
                fld_list = glob.glob('*.maps.fld')
            else:
                return
            
        for fld in fld_list:
            fileList.append(fld.split('.')[0])
        entryfield_value = ""
        if fileList:
            fileList.sort()
            entryfield_value = fileList[0]
    
        ifd = InputFormDescr(title="Run AutoLigand")

        ifd.append({'name':"fileName",
                    'widgetType':Pmw.Group,
                    'container':{'fileName':'w.interior()'},
                    'wcfg':{'tag_text':"Base Name for Receptor Map Files"},
                    })       
        ifd.append({'widgetType':Pmw.ComboBox,
                    'name':'FileBaseName',
                    'parent': 'fileName',
                    'tooltip':'FileBaseName = Just the name part from map files (i.e., FileBaseName.C.map)',
                    'wcfg':{
                            'dropdown':1,
                            'scrolledlist_items':fileList,
                            'entryfield_value':entryfield_value,
                            'selectioncommand':self.selectGrid,
                            'labelpos':'w',}
                    })
        ifd.append({'name':"nGroup",
                    'widgetType':Pmw.Group,
                    'container':{'nGroup':'w.interior()'},
                     'tooltip':'Number of heavy atom for ligand or number of fill points to use.',
                    'wcfg':{'tag_text':"Number of "},
                    })          
        ifd.append({'widgetType':Pmw.RadioSelect,
                    'parent':'nGroup',
                    'name':'type',
                    'defaultValue':'Points',
                    'listtext':["Heavy Atoms","Points"],
                    })        
        ifd.append({'widgetType':Pmw.EntryField,
                    'parent':'nGroup',
                    'name':'number',
                    'wcfg':{
                            'labelpos':'w',
                            'value':'100',
                            'validate':{'validator':'integer'}
                            }
                    })

        ifd.append({'widgetType':Pmw.NoteBook,
                                'name':'autoligandNotebook',
                                'container':{ 'Single Fill':"w.page('Single Fill')",
                                                    'Connected Fill':"w.page('Connected Fill')",
                                                    'Scan':"w.page('Scan')"
                                                  },
                                'componentcfg':[ {'name':'Single Fill', 'cfg':{}},
                                                           {'name':'Connected Fill', 'cfg':{}},
                                                           {'name':'Scan','cfg':{}},
                                                   #        {'raisecommand':self.pageChanged}
                                                         ],
                                'wcfg':{'raisecommand':self.pageChanged},
                                'gridcfg':{'sticky':'we'},
                                })        

        #Single Fill
        ifd.append({'name':"StartLoc",
                    'parent':'Single Fill',
                    'widgetType':Pmw.Group,
                    'container':{'StartLoc':'w.interior()'},
                    'wcfg':{'tag_text':"Start Location"},
                    })        
        ifd.append({'widgetType':ExtendedSliderWidget,
                'name':'gridPointsX',
                'parent':'StartLoc',
                'wcfg':{'label':'X: ',
                        'width':190,
                        'immediate':1,
                        'command':self.changeStartCross,
                        'entrypackcfg':{'side':'left'},},
                })        
        ifd.append({'widgetType':ExtendedSliderWidget,
                'name':'gridPointsY',
                'parent':'StartLoc',
                'wcfg':{'label':'Y: ',
                        'width':190,
                        'immediate':1,
                        'command':self.changeStartCross,
                        'entrypackcfg':{'side':'left'},},
                })        
        ifd.append({'widgetType':ExtendedSliderWidget,
                'name':'gridPointsZ',
                'parent':'StartLoc',
                'wcfg':{'label':'Z: ',                    
                        'width':190,
                        'immediate':1,                        
                        'command':self.changeStartCross,                        
                        'entrypackcfg':{'side':'left'},},
                })        

        ifd.append({'name':"output",
                    'widgetType':Pmw.Group,
                    'parent':'Single Fill',
                    'container':{'output':'w.interior()'},
                    'wcfg':{'tag_text':"Output Options"},
                    })        
        ifd.append({'name':'pdbFile',
                    'tooltip':"""Creates PDB_fill_#Nout1.pdb file where #N is the number of fill points.""",
                    'parent':'output',
                    'widgetType':Tkinter.Checkbutton,
                    'defaultValue': 1,
                    'wcfg':{'text':'Create PDB of the Final Fill',
                            'state':'disabled',
                            
                            },
                    'gridcfg':{'sticky':'w'}
                    })        

        ifd.append({'name':'showProgress',
                    'parent':'output',
                    'tooltip':"""Save intermediate results in a file and open flood player when AutoLigand finishes.""",
                    'widgetType':Tkinter.Checkbutton,
                    'defaultValue': 0,
                    'wcfg':{'text':'Save Intermediate Results for Movie',
                            'variable':Tkinter.IntVar(),
                            },
                    'gridcfg':{'sticky':'w'}
                    })        
        #Connected Fill
        ifd.append({'widgetType':Tkinter.Label,
                'name':'label',
                'parent':'Connected Fill',
                'wcfg':{'text':'Select End Location for Connected Fill'}
                        })   
        ifd.append({'name':"EndLoc",
                    'parent':'Connected Fill',
                    'widgetType':Pmw.Group,
                    'container':{'EndLoc':'w.interior()'},
                    'wcfg':{'tag_text':"End Location"},
                    })        
        ifd.append({'widgetType':ExtendedSliderWidget,
                'name':'gridPointsXEnd',
                'parent':'EndLoc',
                'wcfg':{'label':'X: ',
                        'width':190,
                        'immediate':1,
                        'command':self.changeEndCross,
                        'entrypackcfg':{'side':'left'},},
                })        
        ifd.append({'widgetType':ExtendedSliderWidget,
                'name':'gridPointsYEnd',
                'parent':'EndLoc',
                'wcfg':{'label':'Y: ',
                        'width':190,
                        'immediate':1,
                        'command':self.changeEndCross,
                        'entrypackcfg':{'side':'left'},},
                })        
        ifd.append({'widgetType':ExtendedSliderWidget,
                'name':'gridPointsZEnd',
                'parent':'EndLoc',
                'wcfg':{'label':'Z: ',                    
                        'width':190,
                        'immediate':1,                        
                        'command':self.changeEndCross,                        
                        'entrypackcfg':{'side':'left'},},
                })
                
        #Connected Fill
#        ifd.append({'widgetType':Tkinter.Label,
#                'name':'label',
#                'parent':'Scan',
#                'wcfg':{'text':"""Results from the Scan mode will 
#be available under\n"""+os.getcwd()}
#                        })   
        ifd.append({'widgetType':Pmw.EntryField,
                    'parent':'Scan',                    
                    'name':'nFills',
                    'wcfg':{'label_text':"Number of Fills:",
                            'labelpos':'w',
                            'value':'10',
                            'validate':{'validator':'integer'}
                            }
                    })
                
        def initselect(arg):
            self.selectGrid(entryfield_value)

        self.ifd = ifd
        self.save = self.vf.ICmdCaller.commands.value["Shift_L"]
        self.vf.setICOM(self, modifier="Shift_L", topCommand=0)
        self.vf.setIcomLevel( Atom )
        self.endCross.Set(visible=0)
        val = self.vf.getUserInput(ifd,modal=0, blocking=1,
                                   initFunc=initselect)
        if self.save:
            self.vf.setICOM(self.save, modifier="Shift_L", topCommand=0)
            self.save = None
       # import pdb;pdb.set_trace()
        if val: 
            if not val['FileBaseName'][0]:
                msg = "AutoGrid files are missing.\n"
                msg += "Please generate grid maps. AutoLigand requires input AutoGrid maps."
                tkMessageBox.showerror("Error!", msg)
                return
            selection =  self.ifd.entryByName['autoligandNotebook']['widget'].getcurselection()
            cmdString = [sys.executable, AutoLigandPath]
            cmdString.append('-r')
            cmdString.append(val['FileBaseName'][0])
            self.fileBaseName = val['FileBaseName'][0]
            if val['type'] == 'Points':
                cmdString.append('-p')
                cmdString.append(str(val['number']))
            else:
                cmdString.append('-a')
                cmdString.append(str(val['number']))    
                val['number'] = 6*int(val['number'])                

            if selection in ['Single Fill','Connected Fill']:
                cmdString.append('-x')
                cmdString.append(str(val['gridPointsX']))
                cmdString.append('-y')
                cmdString.append(str(val['gridPointsY']))
                cmdString.append('-z')
                cmdString.append(str(val['gridPointsZ']))
                if selection == 'Connected Fill':
                    cmdString.append('-i')
                    cmdString.append(str(val['gridPointsXEnd']))
                    cmdString.append('-j')
                    cmdString.append(str(val['gridPointsYEnd']))
                    cmdString.append('-k')
                    cmdString.append(str(val['gridPointsZEnd']))
                if os.name == 'nt': #sys.platform == "win32":
                    self.cmdTxt = subprocess.list2cmdline(cmdString)
                else:
                    self.cmdTxt = ' '.join(cmdString)            
                if val['showProgress']:
                    self.cmdTxt += " -m"
                    self.showPlayer = True
                else:
                    self.showPlayer = False
                self.cmd = SysCmdInThread(self.cmdTxt, shell=True)
                self.cmd.start()
                self.checkResults()            
                self.N_of_pts = val['number']
            elif selection == "Scan":
                cmdString.append('-f')
                cmdString.append(str(val['nFills']))
                if os.name == 'nt': #sys.platform == "win32":
                    self.cmdTxt = subprocess.list2cmdline(cmdString)
                else:
                    self.cmdTxt = ' '.join(cmdString)            
                    self.cmdTxt += ' &'#run in a background
                
                subprocess.Popen(self.cmdTxt, shell=True)
                tkMessageBox.showinfo("Success!", """AutoLigand launched successfully.""")                              
        else:
            self.hideGeoms()
            
    def pageChanged(self, name):
        if name == 'Single Fill':
            self.startCross.Set(visible=1)
        elif name == 'Connected Fill':
            self.endCross.Set(visible=1)
        else:
            self.startCross.Set(visible=0)
            self.endCross.Set(visible=0)            
    
    def selectGrid(self, value):
        if not value:
            return
        lines = open(value+".maps.fld").readlines()
        if not lines: return
        for line in lines:
            if line.startswith("#SPACING"):
                spacing = float(line.split()[-1])
            if line.startswith("#NELEMENTS"):
                tmp = line.split()
                dimX = int(tmp[1])
                dimY = int(tmp[2])
                dimZ = int(tmp[3])
            if line.startswith("#CENTER"):
                tmp = line.split()
                centerX = float(tmp[1])
                centerY = float(tmp[2])
                centerZ = float(tmp[3])
        #this variables are the same used in AutoLigand.py
        self.xcent = int(dimX/2)
        self.ycent = int(dimY/2)
        self.zcent = int(dimZ/2)        
        self.centerx = centerX
        self.centery = centerY
        self.centerz = centerZ
        self.spacing = spacing
        c = [centerX,centerY,centerZ]
        xlen = round(spacing*dimX, 4)
        ylen = round(spacing*dimY, 4)
        zlen = round(spacing*dimZ, 4)  
        self.minX = c[0]-xlen*0.5
        self.maxX = c[0]+xlen*0.5
        if self.grids.has_key(value):
            centerX = self.grids[value][0]
            centerY = self.grids[value][1]
            centerZ = self.grids[value][2]

        self.ifd.entryByName['gridPointsX']['widget'].configure(minval=self.minX,
                                                                maxval=self.maxX, )
        self.ifd.entryByName['gridPointsXEnd']['widget'].configure(minval=self.minX,
                                                                maxval=self.maxX, )
        
        self.minY = c[1]-ylen*0.5
        self.maxY = c[1]+ylen*0.5

        self.ifd.entryByName['gridPointsY']['widget'].configure(minval=self.minY,
                                                                maxval=self.maxY, )
        self.ifd.entryByName['gridPointsYEnd']['widget'].configure(minval=self.minY,
                                                                maxval=self.maxY, )
        self.minZ = c[2]-zlen*0.5
        self.maxZ = c[2]+zlen*0.5
        self.ifd.entryByName['gridPointsZ']['widget'].configure(minval=self.minZ,
                                                                maxval=self.maxZ, )
        self.ifd.entryByName['gridPointsZEnd']['widget'].configure(minval=self.minY,
                                                                maxval=self.maxY, )
        
        self.startCross.Set(vertices=((centerX,centerY,centerZ),))        
        self.ifd.entryByName['gridPointsX']['widget'].set(centerX)
        self.ifd.entryByName['gridPointsY']['widget'].set(centerY)
        self.ifd.entryByName['gridPointsZ']['widget'].set(centerZ)
        self.ifd.entryByName['gridPointsXEnd']['widget'].set(centerX+12)
        self.ifd.entryByName['gridPointsYEnd']['widget'].set(centerY+12)
        self.ifd.entryByName['gridPointsZEnd']['widget'].set(centerZ+12)
        self.endCross.Set(vertices=((centerX+12,centerY+12,centerZ+12),))

        pts = [     (self.maxX, self.maxY, self.minZ),
                    (self.minX, self.maxY, self.minZ),
                    (self.minX, self.minY, self.minZ),
                    (self.maxX, self.minY, self.minZ),
                    (self.maxX, self.maxY, self.maxZ),
                    (self.minX, self.maxY, self.maxZ),
                    (self.minX, self.minY, self.maxZ),
                    (self.maxX, self.minY, self.maxZ)
                    ]
         
        self.box.Set(visible=1)
        #self.startCross.Set(visible=1)
        self.box.vertexSet.vertices.array[:] = pts
        self.box.RedoDisplayList()
        self.vf.GUI.VIEWER.Normalize_cb()
        self.vf.GUI.VIEWER.Redraw()
        
    def changeStartCross(self, val):
        centerX = self.ifd.entryByName['gridPointsX']['widget'].get()
        centerY = self.ifd.entryByName['gridPointsY']['widget'].get()
        centerZ = self.ifd.entryByName['gridPointsZ']['widget'].get()
        grid = self.ifd.entryByName['FileBaseName']['widget'].get()
        if grid:
            self.grids[grid] = [centerX,centerY,centerZ]
        self.startCross.Set(visible=1)
        self.startCross.Set(vertices=((centerX,centerY,centerZ),))
        self.vf.GUI.VIEWER.Redraw()   

    def changeEndCross(self, val):
        centerX = self.ifd.entryByName['gridPointsXEnd']['widget'].get()
        centerY = self.ifd.entryByName['gridPointsYEnd']['widget'].get()
        centerZ = self.ifd.entryByName['gridPointsZEnd']['widget'].get()
        grid = self.ifd.entryByName['FileBaseName']['widget'].get()
        #self.endCross.Set(visible=1)
        self.endCross.Set(vertices=((centerX,centerY,centerZ),))
        self.vf.GUI.VIEWER.Redraw()   

    def checkResults(self):
        """Checks the queue for results until we get one"""
        if self.cmd.ok.configure()['state'][-1] == 'normal':
            if self.showPlayer:
                self.openPklData()
            else:
                molName = "FILL_"+str(self.N_of_pts)+"out1"
                if os.path.exists(molName+".pdb"):
                    self.vf.readMolecule(molName+".pdb")
                    self.vf.displaySticksAndBalls(molName, cradius=0.0, sticksBallsLicorice="Sticks and Balls")
                    self.vf.displayLines(molName, negate=True, displayBO=False)
                    self.vf.colorByAtomType(molName, ['sticks', 'balls'])
                self.vf.GUI.ROOT.after(2050, self.hideGeoms)
            return
        self.vf.GUI.ROOT.after(100, self.checkResults)
        
    def hideGeoms(self):
        self.box.Set(visible=0)
        self.startCross.Set(visible=0)
        self.endCross.Set(visible=0)        
        self.spheres.Set(visible=0)
        self.halo.Set(visible=0)
        
    def __call__(self, atom, **kw):
        if not atom:
            return 'ERROR'
        atoms = self.vf.expandNodes(atom)
        if not atoms:
            return 'ERROR'
        atoms = atoms.findType(Atom)
        if not atoms:
            return 'ERROR'
        apply( self.doitWrapper, (atoms,), kw)

    def doit(self, atoms=None):
        if len(atoms)==0: return 
        atom = atoms[0]
        if self.minX < atom.coords[0] > self.maxX:
            return
        if self.minY < atom.coords[1] > self.maxY:
            return
        if self.minZ < atom.coords[2] > self.maxZ:
            return
        
        selection =  self.ifd.entryByName['autoligandNotebook']['widget'].getcurselection()
        if selection ==  "Single Fill":
            self.ifd.entryByName['gridPointsX']['widget'].set(atom.coords[0])
            self.ifd.entryByName['gridPointsY']['widget'].set(atom.coords[1])
            self.ifd.entryByName['gridPointsZ']['widget'].set(atom.coords[2])
    #        self.startCross.Set(visible=1)
            self.startCross.Set(vertices=((atom.coords[0], atom.coords[1], atom.coords[2]),))
        elif selection ==  "Connected Fill":
            self.ifd.entryByName['gridPointsXEnd']['widget'].set(atom.coords[0])
            self.ifd.entryByName['gridPointsYEnd']['widget'].set(atom.coords[1])
            self.ifd.entryByName['gridPointsZEnd']['widget'].set(atom.coords[2])
            self.endCross.Set(vertices=((atom.coords[0], atom.coords[1], atom.coords[2]),))
            
        self.vf.GUI.VIEWER.Redraw()                          
         
        
    def openPklData(self):
        self.vf.AutoLigandMovieCommand(self.fileBaseName+'_flood.pkl')
            
AutoLigandCommandGUI = CommandGUI()
AutoLigandCommandGUI.addMenuCommand('menuRoot', 'Compute', 'Run AutoLigand...', 
                                    cascadeName="AutoLigand")


class OpenMovieCommand(MVCommand):        
    "GUI that opens saved movie of fill file."       
    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        self.floodFile = None


    def guiCallback(self):
        file = self.vf.askFileOpen(types=[('Saved Flood File', '*.pkl'),('All Files', '*')],
                                   title='Open Saved Flood File')        
        if file:
            self.doitWrapper(file)

    def doit(self, file):
        FloodPlayer(self, file)
    
    
OpenMovieCommandGUI = CommandGUI()
OpenMovieCommandGUI.addMenuCommand('menuRoot', 'Compute', 'Open Saved Movie...', 
                                    cascadeName="AutoLigand")

commandList  = [{'name':'AutoLigandCommand','cmd':AutoLigandCommand(),'gui':AutoLigandCommandGUI},
                {'name':'AutoLigandMovieCommand','cmd':OpenMovieCommand(),'gui':OpenMovieCommandGUI}]

def initModule(viewer):
    for _dict in commandList:
        viewer.addCommand(_dict['cmd'],_dict['name'],_dict['gui'])

