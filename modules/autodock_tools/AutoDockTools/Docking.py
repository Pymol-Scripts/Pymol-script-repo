#############################################################################
#
# Author: Ruth HUEY, William  LINDSTROM
#
# Copyright: M. Sanner TSRI 2002
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Docking.py,v 1.41.4.2 2009/09/18 17:57:28 rhuey Exp $
#
# $Id: Docking.py,v 1.41.4.2 2009/09/18 17:57:28 rhuey Exp $
#
#
#
#
#
#
#

"""
This Object is the result of an AutoDock job. 

"""
import os
from string import replace, rfind, find, strip

from MolKit.pdbParser import PdbqParser, PdbqtParser
from MolKit.protein import Residue, ResidueSet




from AutoDockTools.DockingParameters import DockingParameters
from AutoDockTools.DlgParser import DlgParser
from AutoDockTools.Conformation import Conformation, ConformationHandler
from AutoDockTools.Conformation import PopulationHandler
from AutoDockTools.cluster import Clusterer


class Docking:
    """ entity built from an AutoDock docking
        dlo: docking log object built from one docking log file
        parser: docking log parser
        dpo: docking parameter object
        ligMol: ligand molecule
        ligFile: ligand filename

    """

    def __init__(self, parser=None):
        """

        """
        self.defaultParser = 0
        if not parser:
            self.defaultParser = 1
            parser = DlgParser()
        self.parser = parser
        self.dlo_list = []
        self.warnings = []
        self.flex_res = ResidueSet()
        self.ch = None
        self.ph = None   #PopulationHandler
        self.clusterer_dict = {}


    def readDlg(self, dlgFile, ligand=None):
        #only pass dlgFile to init of DockingLogObject 
        #if you want to parse it with DlgParser
        if ligand is not None:
            self.ligMol = ligand
        if self.defaultParser:
            dlo = DockingLogObject(self, dlgFile)
        else:
            dlo = DockingLogObject(self, dlgFile, self.parser)

        if len(dlo.parser.WARNINGS):
            self.warnings.extend(dlo.parser.WARNINGS)

        self.version = dlo.version
        if not self.ch:
            self.ch = ConformationHandler(self.ligMol, 
                                      dlo.dpo['about']['value'])

        oldIndex = len(self.ch.conformations)                              

        #this will just add more conformations
        self.addConformations(dlo.parser)

        if hasattr(dlo.parser, 'population') and len(dlo.parser.population):
            self.ph = PopulationHandler(self.ligMol, 
                                      dlo.dpo['about']['value'])
            for pop in dlo.parser.population_list:
                self.addPopulation(pop, dlo.parser)

        dlo.conformations = self.ch.conformations[oldIndex:]
        # rebuild clusters if there are any
        dlo.hasClusters = 0
        #FIX THIS: always build a clusterer 
        
        cluster_type = "docking"
        if self.version>=4.0:
            cluster_type = "binding"
            
        dlo.clusterer = Clusterer(dlo.conformations, sort=cluster_type)
        #dlo.clusterer = Clusterer(dlo.conformations, sort='docking')
        if dlo.parser.clusterRecord:
            dlo.hasClusters = 1
            #by default sort is on; it just builds argsort list
            #dlo.clusterer = Clusterer(dlo.conformations, sort=0)
            rmstol = dlo.dpo['rmstol']['value']
            #FIX THIS: 
            try:
                dlo.clusterer.rebuild_clusters(dlo.parser.clusterRecord,
                                               rmstol)
                #NB: any dlg has only 1 rmstol value
                dlo.cluster_list = dlo.clusterer.clustering_dict[rmstol]
            except:
                print 'unable to rebuild_clusters '
            rmsref = dlo.dpo['rmsref']['value']
            if len(rmsref):
                dlo.rmsref = rmsref
        self.dlo_list.append(dlo)
        len_dlo_list = len(self.dlo_list)
        if len_dlo_list==1:
            clusterer = dlo.clusterer
        else:
            #each time you add some conformations, need a new clusterer
            clusterer = Clusterer(self.ch.conformations, sort=cluster_type)
            #clusterer = Clusterer(self.ch.conformations, sort='docking')
        self.clusterer = clusterer
        self.clusterer_dict[cluster_type] = clusterer
        #self.clusterer_dict['docking'] = clusterer
        

    def readEntropiaResults(self, dpf, resultsfile):
        # an entropia result involves > 1 resultsfile????
        # parse the results file
        self.parser.parse(resultsfile)
        dlo = DockingLogObject(self)
        dlo.parser = self.parser
        self.dlo_list.append(dlo)
        # parse the dpf file
        dlo.dpo = DockingParameters()
        dlo.dpo.read(dpf)
        dlo.macroStem = dlo.dpo.receptor_stem
        dlo.macroFile = dlo.macroStem + '.pdbqs'
        ####?????####
        self.version = 4
        from MolKit import Read
        # if molecule not found, set ligMol to None 
        try:
            self.ligMol = Read(dlo.dpo['move']['value'])[0]
            # handle the conformations
            if not self.ch:
                #NB all dpo must have the same 'about' and 'move' values
                # and therefore the same ligMol
                self.ch = ConformationHandler(self.ligMol, 
                                      dlo.dpo['about']['value'])
            oldIndex = len(self.ch.conformations)
            self.addConformations(self.parser)
            #dlo.ligMol = self.ligMol
            dlo.conformations = self.ch.conformations[oldIndex:]
        except IOError:
            self.ligMol = None
            dlo.conformations = []


    def set_ligand(self, ligand, dlo=None):
        self.ligMol = ligand
        if dlo is None:
            dlo = self.dlo_list[0]
        if not self.ch:
            #NB all dpo must have the same 'about' and 'move' values
            # and therefore the same ligMol
            self.ch = ConformationHandler(self.ligMol, 
                                  dlo.dpo['about']['value'])
        oldIndex = len(self.ch.conformations)
        #print "oldIndex =", oldIndex
        #IS THIS NECESSARY HERE???
        self.addConformations(self.parser)
        #dlo.conformations = self.ch.conformations
        dlo.conformations = self.ch.conformations[oldIndex:]
        #at the moment this is necessary for making clusterings
        for conf in self.ch.conformations:
            conf.energy = conf.binding_energy


    def readXMLResults(self, resultsfile, build_ligand=True, 
                        dpf=None, ligMol=None):
        # parse the xml results file
        #???assert self.parser.__class__==XMLParser
        self.parser.parse(resultsfile)
        #FORCE VERSION to 4.0
        self.version = self.parser.version = 4.0 ###???###
        dlo = DockingLogObject(self)
        dlo.parser = self.parser
        dlo.filename = dpf
        dlo.dlgFile = resultsfile
        dlo.version = self.parser.version
        dlo.output = None
        dlo.hasClusters = 0
        if dpf is None:
            dpf = self.parser.dpf
        fptr = open(dpf)
        dpflines = fptr.readlines()
        fptr.close()
        dlo.dpo = dlo._buildDpo(dpflines)
        dlo.macroStem = dlo.dpo.receptor_stem
        dlo.macroFile = dlo.macroStem + '.pdbqt'
        if build_ligand:
            ligFile = dlo.dpo['move']['value']
            #try to build the ligand
            try:
                lptr = open(ligFile)
                liglines = lptr.readlines()
                lptr.close()
                ligMol = dlo._buildInputLig(liglines, ligFile)
                ligMol.filename = ligFile
                import os
                ligMol.name = os.path.splitext(ligMol.filename)[0]
                self.set_ligand(ligMol, dlo)
            except:
                print "unable to read ligand file: ", ligFile
                self.ligMol = None
                dlo.conformations = []
        else:
            self.ligMol = ligMol
            if ligMol is not None:
                self.set_ligand(ligMol, dlo)
        self.dlo_list.append(dlo)
        #print "len(dlo_list)=", len(self.dlo_list)


    def addConformations(self, parser):
        self.ch.add(parser.clist, parser.keywords, filename=parser.filename)

    def addPopulation(self, population, parser):
        if not hasattr(self, 'ph'): 
            print 'current docking does not have a population handler'
            return 
        self.ph.add(population, parser.keywords)

    def write_current_conformation(self, filename=""):
        if not hasattr(self, 'ligMol'):
            print "docking has not ligand molecule!"
            return "ERROR"
        liglines = self.ligMol.parser.allLines
        if filename=="":
            ind = str(self.ch.conformations.index(self.ch.current_conf))
            dlgfile= self.dlo_list[0].parser.filename
            dlg_stem = os.path.splitext(os.path.basename(dlgfile))[0]
            parser = self.ligMol.parser
            extension = os.path.splitext(self.ligMol.parser.filename)[1]
            filename = dlg_stem + "_conf" + ind + extension
        fptr = open(filename, 'w')
        ctr = 0
        for l in liglines:
            if l.find("ATOM")!=0 and l.find("HETATM")!=0:
                l += "\n"
                fptr.write(l)
            else:
                crds = self.ligMol.allAtoms[ctr].coords 
                rec = "%s%8.3f%8.3f%8.3f%s\n"%(l[:30],crds[0], crds[1], crds[2],l[54:] ) 
                fptr.write(rec)
                ctr += 1
        fptr.close()
    

    def writeConformation(self, conf=0, filename="", restore=1):
        # if conf==0: write lowest energy conformation
        #this includes setting ligMol to conf and optionally restoring it
        if self.ligMol is None:
            print "Docking has no ligand molecule: unable to write conformation"
            return
        dlo = self.dlo_list[0]
        if conf==0 and hasattr(dlo, 'cluster_list'):
            #CHECK THIS: does dlo_list[0] ALWAYS have a cluster_list??
            conf = dlo.cluster_list[0][0]
        #be sure conf is of class "Conformation"
        if conf.__class__!= Conformation:
            print "conformation must be an instance of Conformation class"
            return
        old_conf = self.ch.current_conf
        self.ch.set_conformation(conf)
        newLines = []
        coords = self.ligMol.allAtoms.coords
        allLines = self.ligMol.parser.allLines
        ext = os.path.splitext(self.ligMol.parser.filename)[-1]
        ctr = 0
        for l in allLines:
            if l[-1]!='\n': l = l + '\n'
            if find(l, 'ROOT')==0 and hasattr(dlo.parser, 'clusterRecord') and dlo.parser.clusterRecord is not None:
                #put the RMS_REF REMARK here
                rms_ref = dlo.parser.clusterRecord[0][0][5]
                nl = "REMARK RMS_REF % 6.4f\n" % (rms_ref)
                newLines.append(nl)
                newLines.append(l)
            elif find(l, 'ATOM')==0 or find(l, 'HETA')==0:
                cc = coords[ctr]
                restLine = l[54:]
                newLines.append(l[:30]+'%8.3f%8.3f%8.3f'%(cc[0],cc[1],cc[2])+l[54:])
                ctr = ctr + 1
            else:
                newLines.append(l)
        if not len(filename):
            filename = self.ligMol.name + "_conf_" + \
                        int(self.ch.conformations.index(conf)) + ext
        fptr = open(filename, 'w')
        for l in newLines:
            fptr.write(l)
        fptr.close()
        if restore and old_conf is not None:
            self.ch.set_conformation(old_conf)
                    



class DockingLogObject:
    """ entity built from parsing a single AutoDock docking
        parser: docking log parser
        dpo: docking parameter object
        docking: Docking this dlo 
        ligMol: ligand molecule
        ligFile: ligand filename

    """

    def __init__(self, docking, filename=None, parser=None):
        self.docking = docking
        self.filename = filename
        if filename:
            if not parser:
                parser = DlgParser()
            self.parser = parser
            parser.parse(filename)
            #for now, keep info about ad3 vs ad4
            self.version = parser.version
            self.dlgFile = filename
            self.dpo = self._buildDpo(parser.dpfLines)
            #if the autodock job did not finish, all models do not exist
            #so shorten the expected list of conformations to actual value
            if len(parser.modelList)<self.dpo['ga_run']['value']:
                self.dpo['ga_run']['value']= len(parser.modelList)
            self.output = parser.histogramlines
            ##???
            self.macroStem = self.dpo.receptor_stem
            self.macroFile = self.macroStem + '.pdbqs'
            if self.parser.version>=4.0:
                self.macroFile = self.macroStem + '.pdbqt'
            ##???
            filename = self.dpo['move']['value']
            if self.docking and not hasattr(self.docking, 'ligMol'):
                ligMol = self._buildInputLig(parser.ligLines, filename=filename)
                if ligMol is not None:
                    ligMol.filename = self.dpo['move']['value']
                    ligMol.name = os.path.splitext(ligMol.filename)[0]
                    self.docking.ligMol = ligMol
                    ligMol.trueLigAtoms = ligMol.allAtoms
                    if hasattr(self, 'flex_res') and self.flex_res:
                        self.docking.flex_res = self.flex_res
                        trueLigAtoms = ligMol.trueLigAtoms
                        for r in self.flex_res:
                            if r.name in ligMol.chains.residues.name:
                                for a in r.atoms:
                                    thisAt = ligMol.allAtoms.get(a.name)
                                    trueLigAtoms = trueLigAtoms - thisAt
                        ligMol.trueLigAtoms = trueLigAtoms 
                    ligMol.lenNonHAtoms = len(ligMol.trueLigAtoms.get(lambda x: x.element!='H'))



    def _buildDpo(self, lines):
        dpo = DockingParameters()
        dpo._parse(lines)
        return dpo


    def _hackLigLineFormat(self, lines):
        for i in range(len(lines)):
            l = lines[i]
            if find(l, 'ATOM')>-1:
                break
        key = lines[i][60:66]
        if len(strip(key)):
            try:
                float(key)
                print 'passed check for extra spaces'
                return lines
            except:
                print 'failed check for extra spaces'
                nl = lines[:i]
                for l in lines[i:]:
                    if find(l,'ATOM')>-1:
                        #this stinks!!!
                        nl.append(l[:54]+l[57:60] + l[61:])
                    else:
                        nl.append(l)
                return nl
        else:
            #can't use this hack
            return lines

            
    def process_flex_res_lines(self, lines):
        # processremove BEGIN_RES, REMARK lines
        # change ROOT/ENDROOT to BRANCH/ENDBRANCH LigROOTATOM# first resAtom#
        #print "in pfrl: len(lines)=", len(lines)
        parser = PdbqtParser()
        pl = []
        flex_res = ResidueSet()
        for i in range(len(lines)):
            l = lines[i]
            #print l
            if l.find('BEGIN_RES')>-1: 
                #print "found BEGIN_RES", len(pl)
                pl = [l]
            elif l.find('END_RES')>-1: 
                #print "found END_RES", len(pl)
                if len(pl):
                    #pl.append(l)
                    parser.allLines = pl[1:-1]
                    mol = parser.parse()[0]
                    mol.chains.residues[-1].torTree = mol.torTree
                    #for res in mol.chains.residues:
                    flex_res += mol.chains.residues
                    #flex_res.extend(mol.chains.residues)
                    #print '2:', len(flex_res)
                    if i!=len(lines)-1:
                        pl = []
            else:
                pl.append(l)
        #print "returning ", flex_res
        return flex_res
                
        
    def _buildInputLig(self, lines, filename):
        if self.parser.version>=4.0:
            parser = PdbqtParser()
            #pdbqParser.allLines = self._hackLigLineFormat(lines)
        else:
            parser = PdbqParser()
        if hasattr(self.parser, 'hasFlexRes') and self.parser.hasFlexRes:
            lines.extend(self.parser.flex_res_lines)
            self.flex_res = self.process_flex_res_lines(self.parser.flex_res_lines)
            #print "flex_res=", self.flex_res
        parser.allLines = lines
        parser.filename = filename
        #filename = pdbqParser.filename = self.dpo['move']['value']
        #the parser does this
        #nameStr = os.path.splitext(filename)[0]
        try:
            mol = parser.parse()[0]
            return mol
        except:
            print "exception in _buildInputLig"

