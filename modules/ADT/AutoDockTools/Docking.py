#############################################################################
#
# Author: Ruth HUEY, William  LINDSTROM
#
# Copyright: M. Sanner TSRI 2002
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Docking.py,v 1.70.2.7 2011/05/05 19:34:48 rhuey Exp $
#
# $Id: Docking.py,v 1.70.2.7 2011/05/05 19:34:48 rhuey Exp $
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
import os, glob, time, sys
from string import replace, rfind, find, strip

from MolKit.pdbParser import PdbqParser, PdbqtParser
from MolKit.protein import Residue, ResidueSet
from MolKit.molecule import AtomSet, MoleculeSet, HydrogenBond, HydrogenBondSet
from MolKit.stringSelector import CompoundStringSelector
from MolKit import Read




from AutoDockTools.DockingParameters import DockingParameters
from AutoDockTools.DlgParser import DlgParser
from AutoDockTools.Conformation import Conformation, ConformationHandler
from AutoDockTools.Conformation import PopulationHandler
from AutoDockTools.cluster import Clusterer
from AutoDockTools.interactiveHistogramGraph import InteractiveHistogramGraph
from mglutil.math.rmsd import RMSDCalculator
from AutoDockTools.InteractionDetector import InteractionDetector



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
        self.css = CompoundStringSelector()


    def addDockings(self, dlo_list, ligand=None):
        for dlo in dlo_list:
            oldIndex = len(self.ch.conformations)                              
            self.addConformations(dlo.parser)
            if len(dlo.parser.WARNINGS):
                self.warnings.extend(dlo.parser.WARNINGS)
            self.dlo_list.append(dlo)
        cluster_type = "binding"
        if self.version<4.0:
            cluster_type = "docking"
        self.clusterer = Clusterer(self.ch.conformations, sort=cluster_type)
        self.clusterer_dict[cluster_type] = self.clusterer


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
        

    def readVSResult(self, vsResultFile, receptor=None, top=None, modelsAs='conformations'):
        ligandSet = Read(vsResultFile, modelsAs=modelsAs)
        if ligandSet is None:
            return "ERROR in reading vs result file", vsResultFile
        assert len(ligandSet)==1
        ligand = ligandSet[0]
        if receptor is not None:
            #??need to buildBondsByDistance here?
            receptorSet= MoleculeSet([receptor]) #??used??
            receptor.hbond_ct = 0
            rec_name = receptor.name
        #for use further on...
        ligand.hbond_ct = 0
        ligand.buildBondsByDistance()
        #@@parser has version information
        parser = ligand.parser
        hb_list = []
        le_hb_list = []
        lc_hb_list = []
        if parser.isV1:
            #build hydrogen bonds  between receptorSet and ligandSet
            hb_list = parser.hb
            ssStr = "~"
        elif parser.isV2:
            for kk in ['AD_LE_hbd', 'AD_LE_hba']: 
                le_hb_list.extend(parser.__dict__[kk])
            for kk in ['AD_LC_hbd', 'AD_LC_hba']: 
                lc_hb_list.extend(parser.__dict__[kk])
            ssStr = "~~"
        else:
            print "unrecognized parser!"
            return 'ERROR'
        if receptor:  
            for hlist in [hb_list, le_hb_list, lc_hb_list]:
                for hb_str in hlist:
                    #VALIDITY CHECK
                    ok = False
                    #['xJ1_xtal:B:GLY16:N,HN~ZINC02025973_vs:d:<0>:O3\n','xJ1_xtal:B:LEU63:N,HN~ZINC02025973_vs:d:<0>:O5\n']    
                    # v2
                    #d:<0>:N1~~B:ASN83:N,d:<0>:O3~~B:LYS20:NZ
                    str_1, str_2 = hb_str.split(ssStr)
                    if not len(str_1) or not len(str_2):
                        print "skipping improper hb_str ", hb_str
                        continue
                    ligStr = ligand.chains[0].id #@@FIX THIS
                    if parser.isV1:
                        ligStr = ligand.name
                    ok = True
                    if ligStr in str_1 and rec_name in str_2:
                        lig_str = str_1
                        rec_str = str_2
                    elif ligStr in str_2 and rec_name in str_1:
                        lig_str = str_2
                        rec_str = str_1
                    else:
                        print "ligand name: ", ligand.name, " not found in hb_str ", hb_str
                        ok = False                    
                    if ok:
                        #rec_str xJ1_xtal:B:GLY16:N,HN
                        # try to retrieve receptor atoms
                        rec_ats = self.css.select(receptorSet, strip(rec_str))
                        #rec_ats = self.vf.expandNodes(strip(rec_str))
                        if not len(rec_ats):
                            print "skipping hb_str because ",  rec_str, " did not match atoms in the receptor"
                        #(<AtomSet instance> holding 2 Atom, "xJ1_xtal:B:ARG57:NE,HE", '')
                        rec_ats = rec_ats[0]
                        #add ligand name to lig_str
                        lig_ats = self.css.select(ligandSet, strip(lig_str))[0]
                        #(<AtomSet instance> holding 1 Atom, "ZINC02025973_vs_le:d:<0>:O3", '')
                        #lig_ats = self.vf.expandNodes(strip(lig_str))
                        if not len(lig_ats):
                            print "skipping hb_str because ",  lig_str, " did not match atoms in the ligand"
                        #(<AtomSet instance> holding 1 Atom, "ZINC02026663_vs:d:<0>:O3", '')
                        #lig_ats = lig_ats[0]
                        donor_ats = rec_ats
                        accAt = lig_ats[0]
                        if "," in lig_str: 
                            donor_ats = lig_ats
                            accAt = rec_ats[0]
                        if not hasattr(accAt, 'hbonds'):
                            accAt.hbonds = HydrogenBondSet()
                        #?? or ??
                        donAt = donor_ats[0]
                        if not hasattr(donAt, 'hbonds'):
                            donAt.hbonds = HydrogenBondSet()
                        hAt = None
                        if len(donor_ats)==2:
                            hAt = donor_ats[1]
                            if not hasattr(hAt, 'hbonds'):
                                hAt.hbonds = HydrogenBondSet()
                        hb = HydrogenBond(donAt, accAt, hAt)
                        #add them to these atoms
                        donAt.hbonds.append(hb)
                        accAt.hbonds.append(hb)
                        if hAt is not None:
                            if not hasattr(hAt, 'hbonds'):
                                hAt.hbonds = HydrogenBondSet()
                            hAt.hbonds.append(hb)
                        receptor.hbond_ct += 1
                        ligand.hbond_ct += 1
                        #whew!
            #cc
            # parser.lig_close_ats
            #['ZINC02025973_vs:d:<0>:N1,C13,C11,C6,C9,C7,O1,C12,C8,C1,O3,C5,O4,O5']
            if parser.isV1: 
                if len(parser.lig_close_ats) and ligand.name in parser.lig_close_ats[0]:
                    #self.css.select returns atomset, possibly empty PLUS string
                    #(<AtomSet instance> holding 15 Atom, "ZINC00027719_vs:d:<0>:O3,C5,C3...", '')
                    close_ats, close_ats_str = self.css.select(ligandSet, strip(parser.lig_close_ats[0]))
                    if not len(close_ats):
                        print "skipping ligand close ats because ",  parser.lig_close_ats[0], " did not match atoms in the ligand"
                    #ligand.close_ats = self.vf.expandNodes(strip(parser.lig_close_ats[0]))
                    ligand.close_ats = close_ats
                # parser.macro_close_ats
                rec_close_ats = AtomSet()
                if hasattr(parser, 'macro_close_ats') and len(parser.macro_close_ats) and receptor.name in parser.macro_close_ats[0]:
                    for macro_close_ats_str in parser.macro_close_ats:
                        cc_ats = self.css.select(receptorSet, strip(macro_close_ats_str))[0]
                        #cc_ats = self.vf.expandNodes(strip(macro_close_ats_str))
                        if len(cc_ats):
                            rec_close_ats.extend(cc_ats)
                    receptor.close_ats = rec_close_ats.uniq()
                #pi_pi and cation_pi
                ligand.pi_pi = parser.pi_pi
                for pi_str in parser.pi_pi:  
                    ok = False
                    str_1, str_2 = pi_str.split("~~")
                    #contrived example:
                    #pi_st1  'xJ1_xtal:B:TRP42:CD2,CE2,CZ3,CH2,CZ3,CE3,NE1,CD1,CG
                    #pi_st2  'ZINC02026663_vs:d:<0>:C7,C4,C2,C5,C8,N2,N1,C3,C9,C6'
                    if not len(str_1) or not len(str_2):
                        print "skipping improper pi_str ", pi_str
                        continue
                    if ligand.name in str_1 and rec_name in str_2:
                        lig_str = str_1
                        rec_str = str_2
                        ok = True
                    elif ligand.name in str_2 and rec_name in str_1:
                        lig_str = str_2
                        rec_str = str_1
                        ok = True
                    else:
                        print "ligand name: ", ligand.name, " not found in hb_str ", hb_str
                    if ok:
                        #rec_str xJ1_xtal:B:GLY16:N,HN
                        # try to retrieve atoms
                        lig_ats = None
                        rec_ats = None
                        rec_ats = self.css.select(receptorSet, strip(rec_str))
                        #rec_ats = self.vf.expandNodes(strip(rec_str))
                        if not len(rec_ats):
                            print "skipping ", pi_str,  " because ",  rec_str, " did not match atoms in the receptor"
                            continue
                        rec_ats = rec_ats[0]
                        lig_ats = self.css.select(ligandSet, strip(lig_str))
                        #lig_ats = self.vf.expandNodes(strip(lig_str))
                        if not len(lig_ats):
                            print "skipping ", pi_str,  " because ",  lig_str, " did not match atoms in the ligand"
                            continue
                        lig_ats = lig_ats[0]
                        if rec_ats and lig_ats:
                            ligand.pi_pi_ats = lig_ats
                            receptor.pi_pi_ats = rec_ats
                    #ligand.pi_cation = parser.pi_cation
                    for pi_str in parser.pi_cation:  
                        #'xJ1_xtal:B:LYS55:NZ~~ZINC02026663_vs:d:<0>:C7,C4,C2,C5,C8,N2,N1,C3,C9,C6'
                        ok = False
                        str_1, str_2 = pi_str.split("~~")
                        #pi_st1  'xJ1_xtal:B:LYS55:NZ'
                        #pi_st2  'ZINC02026663_vs:d:<0>:C7,C4,C2,C5,C8,N2,N1,C3,C9,C6'
                        if not len(str_1) or not len(str_2):
                            print "skipping improper pi_str ", pi_str
                            continue
                        if ligand.name in str_1 and rec_name in str_2:
                            lig_str = str_1
                            rec_str = str_2
                            ok = True
                        elif ligand.name in str_2 and rec_name in str_1:
                            lig_str = str_2
                            rec_str = str_1
                            ok = True
                        else:
                            print "ligand name: ", ligand.name, " not found in hb_str ", hb_str
                        if ok:
                            # try to retrieve receptor atoms
                            #rec_str xJ1_xtal:B:GLY16:N,HN
                            pi_ats = None
                            cation_at = None
                            rec_ats = self.css.select(receptorSet, strip(rec_str))
                            #rec_ats = self.vf.expandNodes(strip(rec_str))
                            #if not len(rec_ats):
                                #print "skipping ", pi_str,  " because ",  rec_str, " did not match atoms in the receptor"
                            rec_ats = rec_ats[0]
                            if len(rec_ats)>1:
                                receptor.pi_ats = rec_ats 
                                #print "set pi_ats to rec_ats", rec_ats.full_name()
                            elif len(rec_ats)==1:
                                receptor.cation_at = rec_ats[0]
                                #print "set cation_at to rec_ats[0]", rec_ats[0].full_name()
                            lig_ats = self.css.select(ligandSet, strip(lig_str))
                            #lig_ats = self.vf.expandNodes(strip(lig_str))
                            if not len(lig_ats):
                                print "skipping ", pi_str,  " because ",  lig_str, " did not match atoms in the ligand"
                            lig_ats = lig_ats[0] #because css.select returns a tuple
                            if len(lig_ats)==1: 
                                if cation_at is not None:
                                    #print "skipping pi_str because ",  lig_str, " only single atom specified for lig + for rec"
                                    continue
                                else:
                                    cation_at = lig_ats[0]
                                    center1 = cation_at.coords
                                    #print "set cation_at to lig_ats[0]", lig_ats[0].full_name()
                                    ligand.cation_at = cation_at
                            elif pi_ats is not None:
                                print "skipping pi_str because ",  lig_str, " multiple atoms specified for both lig + rec"
                                continue
                            else:
                                ligand.pi_ats = lig_ats
                                #print "set pi_ats to lig_ats", lig_ats.full_name()
            if top is not None:
                #CLUSTERING HISTOGRAM
                if parser.isV1:
                    ligand.cl_lines = parser.cl_lines
                    current_cl = None
                    summary_list = parser.vs_summary_line.split()
                    energy_used = summary_list[1]
                    rms = summary_list[2]
                    num_runs = summary_list[3]
                    #'USER  AD> ligand efficiency  -0.3112\n'
                    #'USER  AD> rmsd     LE   clu:size e_range\n'
                    #ligand_efficiency = ligand.cl_lines[0].split()[-1].strip()
                    #'USER  AD> 0.000  -7.470  * ** 39 0.60 *\n'
                    tstr = ligand.name + ': rms=' + rms + ' ' + num_runs + ' runs'
                    top.title(tstr)
                    #NOW build the histogram
                    xlabel = 'BINDING ENERGY'
                    #dataList=[[-11.23, 8], [-10.380000000000001, 2]]
                    #reverseList = [[0, 1, 2, 3, 4, 5, 6, 7], [8, 9]]
                    dataList=[]
                    reverseList=[]
                    ctr = 2
                    best_energy = False
                    largest_cluster = False
                    other = False
                    num_confs = 0
                    if not len(ligand.cl_lines) or len(ligand.cl_lines)<3:
                        return
                    cl_ctr = 0
                    for cl_line in ligand.cl_lines[2:]:
                        # SKIP FIRST TWO LINES:
                        #'USER  AD> ligand efficiency  -0.5673\n'
                        #'USER  AD> rmsd     LE   clu:size e_range\n'
                        #if ctr<2: continue
                        ll = cl_line.split(',')
                        len_ll = len(ll)
                        # append tuple(lowest_energy, number_in_cluster)
                        #['0.000', '-7.470', '39', '0.60', 'faah8621_ZINC02026663_xJ1_xtal_00.dlg', '30', '1', '1', '1\n']
                        #FIND THE RIGHT LINE FIRST:
                        enrg = round(float(ll[1]),3)
                        num_confs = int(ll[2])
                        new_item = [enrg, num_confs]
                        dataList.append(new_item)
                        if num_confs>1:
                            reverseList.append(range(ctr, ctr+num_confs+1))
                        if ll[-3]=='1':
                            #FOUND IT!!
                            # ONLY interested in best energy and largest cluster
                            if ll[-2].strip()=="1":
                                best_energy = 1
                            elif ll[-1].strip()=="1":
                                largest_cluster = 1
                            else:
                                other = cl_ctr
                                other_energy = enrg
                            #print "dataList=", dataList
                            # [[-8.51, 200],]
                            #[range(0,200)], 
                        ctr += 1
                        cl_ctr += 1
                    ligand.clustNB = InteractiveHistogramGraph(ligand.name,
                        master=top, nodeList = dataList, reverseIndex=reverseList,
                        label_text=ligand.name + ':' + rms + ' rms', xlabel_text=xlabel, 
                        ylabel_text='#\nC\nO\nN\nF\nO\nR\nM\nA\nT\nI\nO\nN\nS')
                    #ligand.clustNB.draw.bind('<Button-1>', CallBackFunction(self.describeClust, ligand),'+')
                    #COLOR ONE BAR RED
                    #IF LOWEST ENERGY:
                    le = +100
                    le_index = None
                    if best_energy:
                        for i,j in ligand.clustNB.geoms.items():
                            if j.point[0]<le:
                                le = j.point[0]
                                le_index = i
                        ligand.clustNB.draw.itemconfig((le_index,), fill='red')
                    max_height = -1
                    max_index = None
                    if largest_cluster:
                        for i,j in ligand.clustNB.geoms.items():
                            if j.height>max_height:
                                max_height = j.height
                                max_index = i
                        ligand.clustNB.draw.itemconfig((max_index,), fill='red')
                    if other:
                        for i,j in ligand.clustNB.geoms.items():
                            if j.point[0]==other_energy:
                                other_index = i
                        ligand.clustNB.draw.itemconfig((other_index,), fill='red')
        self.ligMol = ligand
        return ligandSet

        

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
        num_prev = len(self.ch.conformations)
        self.ch.add(parser.clist, parser.keywords, filename=parser.filename)
        for c in self.ch.conformations[num_prev:]:
            c.filename = os.path.basename(parser.filename)

    def addPopulation(self, population, parser):
        if not hasattr(self, 'ph'): 
            print 'current docking does not have a population handler'
            return 
        self.ph.add(population, parser.keywords)

    def write_current_conformation(self, filename="", infoStr=None, rms=-1, ncl_to_write=-1, summary_only=False):
        #infoStr: if not None include clustering info
        if not hasattr(self, 'ligMol'):
            print "docking has no ligand molecule!"
            return "ERROR"
         #rms: if -1, use smallest key
        if ncl_to_write>0 and rms==-1:
            rms_list = self.clusterer.clustering_dict.keys()
            if not len(rms_list):
                print "docking has no clustering!"
                return "ERROR"
            rms_list.sort()
            rms = rms[0]
        liglines = self.ligMol.parser.allLines
        ind = str(self.ch.conformations.index(self.ch.current_conf))
        if filename=="":
            dlgfile= self.dlo_list[0].parser.filename
            dlg_stem = os.path.splitext(os.path.basename(dlgfile))[0]
            parser = self.ligMol.parser
            extension = os.path.splitext(self.ligMol.parser.filename)[1]
            filename = dlg_stem + "_conf" + ind + extension
        fptr = open(filename, 'w')
        ctr = 0
        if infoStr and hasattr(self, 'clusterer') and len(self.clusterer.clustering_dict):
            #try:
            #    fptr.write(self.clusterer.getInfoStr(infoStr, ind, ncl_to_write, rms))
            #except:
            #    print "error writing clustering information"
            cl_lengths = map(len,self.clusterer.clustering_dict[rms])
            max_ind = cl_lengths.index(max(cl_lengths))
            if ncl_to_write<0:
                ncl_to_write=max_ind+5
            best = self.clusterer.clustering_dict[rms][0][0]
            #@@@
            self.ch.set_conformation(best)
            num_confs = len(self.ch.conformations)
            best_ind = str(self.ch.conformations.index(best))
            # infoStr should be '#'
            if infoStr is None:
                infoStr="#"
            #?IS THIS EVER CALLED??
            sss, dr_str = self.clusterer.getInfoStr(comment=infoStr, ind=best_ind, rms=rms, ncl_to_write=ncl_to_write, include_dlgfilename_run=True)
            sss_list = sss.split('\n')
            dr_str_list = dr_str.split('\n')
            for i in range(len(sss_list)-1):  #starts with '#binding' so numbers are 1-based
                # print clustersize
                if i==0:
                    fptr.write( "%s %4.2f %d \n"%(sss_list[i],rms, num_confs))
                elif i==1:
                    #print energy_range
                    cl = self.clusterer.clustering_dict[2.][0]
                    b_cl = cl[0]
                    w_cl = cl[-1]
                    len_cl = len(cl)
                    #print sss_list[i], "%s %4.2f\n"%(sss_list[i],rms_tolerance))
                    #print " %d %4.2f *"%(len_cl, w_cl.binding_energy - b_cl.binding_energy)
                    fptr.write("%s %d %4.2f %s *\n"%(sss_list[i], len_cl, w_cl.binding_energy - b_cl.binding_energy, dr_str_list[i]))
                elif i == max_ind:
                    m_cl = self.clusterer.clustering_dict[2.][max_ind]
                    mb_cl = m_cl[0]
                    mw_cl = m_cl[-1]
                    #print sss_list[i] + " **",
                    len_cl = len(m_cl)
                    ####print "%4.2f %d"%(mw_cl.binding_energy - mb_cl.binding_energy, len(m_cl) )
                    #print " %d %4.2f *"%(len_cl, mw_cl.binding_energy - mb_cl.binding_energy)
                    fptr.write( "%s %d %4.2f %s *\n"%(sss_list[i], len_cl, mw_cl.binding_energy - mb_cl.binding_energy, dr_str_list[i]))
                else:
                    #print  sss_list[i], len(d.clusterer.clustering_dict[2.0][i])
                    fptr.write("%s %d %s\n"%( sss_list[i], len(self.clusterer.clustering_dict[2.0][i], dr_str_list[i])))
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


from numpy import array, zeros, sqrt  
from glob import glob  
import os
from sys import argv
from mglutil.math.rmsd import RMSDCalculator
from bhtree import bhtreelib
from datetime import date
import time


class FoxResultProcessor:
    """ Create pdbqt-plus file from a directory containing AutoDock results for a single ligand to a single receptor

    """
    atype_list_complete = False
    atype_list = []

    # X-Y hydrogen-bond distance cutoff -X-H...Y-
    HB_CUTOFF=2.75

    # vdw tolerance (account for AutoGrid map smoothing)
    VDW_TOL = .5

    # default RMSD tolerance
    #RMSTOL = 2.0

    DEBUG = True

    rec_acc,rec_don = None, None 

    vdw_radii = { 'H': 1.00, # from AD4_parameter.dat Rii/2 values
                  'HD': 1.00,
                  'HS': 1.00,
                  'C': 2.00,
                  'A': 2.00,
                  'N': 1.75,
                  'NA': 1.75,
                  'NS': 1.75,
                  'OA': 1.60,
                  'OS': 1.60,
                  'F': 1.54,
                  'Mg': 0.65,
                  'MG': 0.65,
                  'P': 2.10,
                  'SA': 2.00,
                  'S': 2.00,
                  'Cl': 2.04,
                  'CL': 2.04,
                  'Ca': 0.99,
                  'CA': 0.99,
                  'Mn': 0.65,
                  'MN': 0.65,
                  'Fe': 0.65,
                  'FE': 0.65,
                  'Zn': 0.74,
                  'ZN': 0.74,
                  'X': 2      } # default vdW for unknown atom


    def __init__(self, receptor_filename=None, RMSTOL=2.,verbose=False):
        self.receptor_filename=receptor_filename
        if receptor_filename==None:
            print "receptor_filename must be specified"
            return "ERROR"
        self.RMSTOL = RMSTOL
        # rec_crds keys are 'text', 'coord', 'atype'
        self.rec_crds = self.get_coords(receptor_filename)
        self.rec_name = os.path.basename(receptor_filename)
        self.rec_mol_name = os.path.splitext(self.rec_name)[0]
        self.verbose = verbose
        if self.verbose:
            print "__init__: len(self.rec_crds['coord'])=", len(self.rec_crds['coord'])


    def get_lines(self, filename):
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        return lines


    def Dist(self, f, s, sq=True):  
        """works with PDB lines"""
        if sq:
            return sqrt((float(f[30:38])-float(s[30:38]))**2 + (float(f[38:46])-float(s[38:46]))**2 +  (float(f[46:54])-float(s[46:54]))**2)
        else:
            return (float(f[30:38])-float(s[30:38]))**2 + (float(f[38:46])-float(s[38:46]))**2 +  (float(f[46:54])-float(s[46:54]))**2


    def getPoses(self, dlg_list, include_hydrogens=False):
        """
        INPUT :parse a list of dlg
        return a list of poses contaning coordinates as text (to be written plantiff) and float array (to be used in calculations)
        poses_list = [
             0: #       text_pose   => original ASCII text of the pose (to be used to write the PDBQT+)
             1: #       pose_coords => numpy float array (to be used to calculate RMSD faster)
             2: #       energy      => float of the energy
        The function initializes also the atom types list (used to calculate vdW interactions later)

        TODO : this function could check the ligand name consistency.
        """
        accepted_kw = [ "ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF", "REMARK" ]
        #global atype_list_complete, atype_list
        ligand_name = None
        poses_list = []
        problematic = []
        for d in dlg_list:
            inside = False
            for l in self.get_lines(d):
                if l[0:7] == "DOCKED:":
                    inside = True
                    l = l[8:]
                    if "MODEL" in l:
                        # open pose
                        text_pose = []
                        coord = []
                        text_pose.append(l)
                    elif "ENDMDL" in l:
                        coord = array( coord, 'f')
                        poses_list.append( { "text" : text_pose, "coord" : coord, "energy" : e } )
                        self.atype_list_complete = True
                        text_pose.append(l)
                    elif l.startswith("ATOM") or l.startswith("HETATM"):
                        type = l.rsplit()[-1]
                        # for the float array calculation HD are excluded 
                        # by default (RMSD calculation *must not* consider HDs!)
                        if (not type == "HD") or include_hydrogens : 
                            #coord.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                            coord.append(map(float, [l[30:38],l[38:46],l[46:54]]))
                            if not self.atype_list_complete: # FUTURE: to be used for the AD David's method
                                self.atype_list.append( type ) 
                            #try:
                            #    coord.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                            #    if not atype_list_complete: # FUTURE: to be used for the AD David's method
                            #        atype_list.append( type ) 
                            #except:
                            #    # problems in parsing these coordinates
                            #    print "WARNING! error in parsing a coords"
                            #    problematic.append(d)
                            #    break 
                        text_pose.append(l)
                    elif "USER    Estimated Free Energy of Binding    =" in l:
                        e = l.split("=")[1]
                        e = float(e.split("k")[0])
                    elif l.split(None, 1)[0] in accepted_kw:
                        text_pose.append(l)
                elif l.startswith("DPF> move") and not ligand_name:
                    ligand_name = l.split("DPF> move ")[1].split(".pdbqt")[0]
                    
                if l.startswith("DPF >") and inside: 
                    # stop reading a dlg file after the latest pose (average lines skipped: ~30%)
                    break
        return poses_list, problematic, ligand_name # poses_list format: [ [text_pose, pose_coords, energy], ...]


    def get_coords(self, filename=None, list=None, include_hydrogens=True):
        """ 
        designed to extract receptor coordinates from a PDB(QT)
        but it can be used with lists too.

        designed mainly to be used for extracting receptor coordinates

        return { 'text':text_pdb_atom_entries, numeric_array_coords, atype]
               { 'text':[atoms], 'coord': array(coord, 'f'),  'atype':[atomtypes]}
        """
        coord = []
        atoms = []
        atype = []
        if not filename and not list:
            return False
        if filename:
            source = self.get_lines(filename)
        if list:
            source = list
        for l in source:
            if l.startswith("ATOM") or l.startswith("HETATM"):
            #if (l[0:4] == "ATOM") or (l[0:6] == "HETATM"):
                at = l.rsplit(None, 1)[1]
                if not at == "HD" or include_hydrogens: # by default, HD are included
                    coord.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                    atoms.append(l)
                    atype.append(at) 
        return { 'text' : atoms, 'coord' : array( coord, 'f'), 'atype': atype }


    def path_to_list(self, path):
        return glob(os.path.join(path, "*.dlg"))


    def get_hb_interactions(self, ligand, receptor):
        """
        input   : pdb lines (ligand, receptor)
        output  : list of acceptor/donor pairs (lig:rec) <= PDBQT+ format
        notes   : distance only is used to characterize HB (consistant with AG maps)
        """
        #global rec_acc, rec_don
        def findHbAccepDon(list = None, filename = None):
            acceptor_types = ['OA', 'NA', 'SA']
            donor_types = ['N', 'O', 'OA', 'NA']
            acceptors = []
            donors = []
            h = []
            dcandidate = []
            H_COV_BOND = 1.1  

            if not list and not filename:
                return
            if list:
                source = list
            if filename:
                source = self.get_lines(filename)
            for l in source:
                if l.startswith("ATOM") or l.startswith("HETATM"):
                    atype=l.split()[-1]
                    if atype in acceptor_types:
                        if not l in acceptors:
                            acceptors.append(l)
                    if atype in donor_types:
                        if not l in dcandidate:
                            dcandidate.append(l)
                    if atype == 'HD':
                        if not l in h:
                            h.append(l)
            H_COV_BOND  = H_COV_BOND ** 2  
            for a in dcandidate:
                for hx in h:
                    if self.Dist(a, hx, sq= False) <= H_COV_BOND:
                        donors.append([a,hx])
            return acceptors, donors # donors are in the form of [ [heavy_atom, hydrogen], [... , ...], ...] 
            
        # use the close contact atoms to calculate the HB (faster)
        lig_acc,lig_don = findHbAccepDon(list = ligand) # [acceptors, donors]
        if not self.rec_acc or not self.rec_don:
            self.rec_acc,self.rec_don = findHbAccepDon(list = receptor)

        # Hydrogen bond identification
        hb_pair = []
        hb_lig_acc = []
        hb_lig_don = []

        # ligand hb acceptor atoms
        for a in lig_acc:
            for r in self.rec_don:
                if self.Dist(a, r[1],sq = False) <= self.HB_CUTOFF**2:
                    if not [a,r] in hb_pair:
                        hb_pair.append([a,r[1]])
                        hb_lig_acc.append([a,r[1]])

        # ligand hb donor atoms
        for a in lig_don:
            for r in self.rec_acc:
                if self.Dist(a[1], r,sq = False) <= self.HB_CUTOFF**2:
                    if not [a,r] in hb_pair:
                        hb_pair.append([a[1],r])
                        hb_lig_don.append([a[1],r])
        return { 'acceptors' : hb_lig_acc, 'donors': hb_lig_don}


    def get_contact_atoms(self, ligand_pose, receptor):
        """
            Input   : ligand pose, receptor
            Output  : list of receptor contact atoms
        """
        #print ligand_pose
        contact_atoms = []
        hb_atoms = []

        cutoff = 7. # initial cutoff for bhtree

        bht = bhtreelib.BHtree( receptor['coord'], None, 10)
        indices = zeros( (len(receptor['coord']),) ).astype('i')
        dist    = zeros( (len(receptor['coord']),) ).astype('f')

        for i in range( len(ligand_pose["coord"])):
            try:
                l_vdw = self.vdw_radii[ self.atype_list[i] ] 
            except:
                l_vdw = self.vdw_radii['X'] 

            nb = bht.closePointsDist(tuple(ligand_pose["coord"][i]), cutoff, indices, dist)
            for j in range(nb):
                rec_index = indices[j]
                rec_atom = receptor['text'][rec_index]
                d = dist[j]
                try:
                    r_vdw = self.vdw_radii[receptor['atype'][rec_index]]
                except:
                    r_vdw = self.vdw_radii['X']
                if d <= ( r_vdw + l_vdw + self.VDW_TOL ): 
                    if not rec_atom in contact_atoms:
                        contact_atoms.append(rec_atom)
                    # TODO potentially interesting to put HB check here # to speed up...?
        return contact_atoms 


    def lig_eff(self, energy, heavy_atoms):
        """ return ligand efficiency calculated on energy/heavyatoms """
        return float(energy)/float(heavy_atoms)
            

    def percent(self, value, total):
        return (float(value)/float(total))*100


    def pmv_atom_stripper(self, mol_name, atom, contract = False ):
        # "ATOM    455  N   GLY A  48      -2.978   5.488   6.818  1.00 11.64    -0.351 N" (xJ1_xtal)
        #
        #       |
        #      \./
        #       '  
        #  xJ1_xtal:A:GLY48:N
        """
        transform pdb entries to PMV syntax
        """
        chain = atom[21].strip()
        res_name = atom[16:21].strip()
        res_num = atom[22:26].strip()
        atom = atom[12:16].strip()
        if not contract:
            return "%s:%s:%s%s:%s" % (mol_name, chain, res_name, res_num, atom)
        else:
            return "%s:%s%s:%s" % (chain, res_name, res_num, atom)


    def FastReclustering(self, poses):
        """
        input   : poses
        output  : poses clustered (as list) => [ [x,x,x], [y,y], [ z,z,z,z,z,z,z,z], ... ]

        """
        clusters = []
        while poses:
            rmsdcalc = RMSDCalculator(poses[0]["coord"])
            current_cluster = [poses[0]]
            del poses[0]
            removable = []
            for p in range(0, len(poses)):
                candidate = poses[p]
                if rmsdcalc.computeRMSD(candidate["coord"]) <= self.RMSTOL:
                    current_cluster.append(candidate)
                    removable.append(candidate)
            for i in removable:
                del poses[poses.index(i)]
            clusters.append(current_cluster)
        return clusters


    def extractResultsPoses(self, clusters):
        # TODO this can be improved
        LC_pose = ""
        LC_size = -1
        LC_energy = 10000.
        LC_index = -1

        LE_pose = ""
        LE_size = -1
        LE_energy = 10000.
        LE_index = -1

        histogram = []

        # identify LE and LC
        c = 0
        for pop in clusters:
            pop.sort(key=lambda x: x['energy'])

            curr_cluster_size = len(pop)
            curr_cluster_energy = pop[0]['energy']

            # populate histogram
            histogram.append([curr_cluster_energy, curr_cluster_size])

            # identify largest cluster # 
            if (curr_cluster_size > LC_size) or (curr_cluster_size == LC_size and curr_cluster_energy < LC_energy ):
                LC_size = curr_cluster_size
                LC_pose = pop[0]
                LC_energy = curr_cluster_energy
                LC_index = c

            if curr_cluster_energy < LE_energy:
                LE_energy = curr_cluster_energy
                LE_pose = pop[0]
                LE_size = curr_cluster_size
                LE_index = c
            c += 1

        # collate results
        results = [ [LE_energy, LE_size, LE_pose]]
        histogram[LE_index].append("**")

        if not LE_pose == LC_pose:
           results.append( [LC_energy, LC_size, LC_pose] )
           histogram[LC_index].append("*")

        # sort the histogram
        histogram.sort(key = lambda x: x[0])
        return results, histogram


    def process(self, receptor_filename=None, directory=None, outputfilename=None, 
                     rms_tolerance=None, dlg_list=[]):
        """
CODE TO REPLACE THESE LINES USED FOR TESTING AND DEVELOPMENT:
    from AutoDockTools.Docking import FoxResultProcessor
    frp = FoxResultProcessor()
    import glob
    fl = glob.glob("xJ1_xtal.pdbqt")
    crds = frp.get_coords(filename="xJ1_xtal.pdbqt")
    dlg_list = glob.glob("faah8621_ZINC02025973_xJ1_xtal/*dlg")
    receptor_filename = "xJ1_xtal.pdbqt"
    poses, b, ligand_name = frp.getPoses(dlg_list)
    total_runs = len(poses)
    # reclustering
    print " => RECLUSTERING ", 
    clusters = frp.FastReclustering(poses)
    results, histogram = frp.extractResultsPoses(clusters)
    pdbqt = frp.createPDBQTplus(results, total_runs, histogram, receptor_filename, ligand_name, dlg_list = dlg_list)

        """
        if directory is None and dlg_list is None:
            print "directory OR dlg_list must be specified"
            return
        if self.verbose:
            print "in process, directory=", directory
            print " len(self.rec_crds)=", len(self.rec_crds), 
            print " len(self.rec_crds['coord']=", len(self.rec_crds['coord'])
        #fl = glob(receptor_filename)
        #assert len(fl)==1
        #self.crds = frp.get_coords(filename=fl[0])
        if directory and not dlg_list:
            dlg_list = glob(os.path.join(directory,"*.dlg"))
        poses, b, ligand_name = self.getPoses(dlg_list)
        total_runs = len(poses)
        clusters = self.FastReclustering(poses)
        results, histogram = self.extractResultsPoses(clusters)
        pdbqt = self.createPDBQTplus(results, total_runs, histogram, receptor_filename, ligand_name, dlg_list = dlg_list, clusters=clusters)
        pdbqt_filename = ligand_name +"_VS.pdbqt"
        if outputfilename is not None:
            pdbqt_filename = outputfilename
        fptr = open(pdbqt_filename, 'w')
        fptr.writelines(pdbqt)
        fptr.close()
        if self.verbose: print "closed ", pdbqt_filename


    def createPDBQTplus(self, results, total_runs, histogram, receptor_filename, lig_mol_name , dlg_list, max_dlg_count = None, clusters=[] ):
        """
            input   :   results list in the form [ 
                            [LE_energy, LE_size, LE_pose], 
                            [LC_energy, LC_size, LC_pose] ]
                        total_runs  = integer
                        histogram   = histogram list 
                        receptor_filename (this will be handled by Fox)
                        lig_mol_name      (this should be handled by Fox)

            output  :   text buffer of the PDBQT+

        """
        # TODO this part will be handled by Fox
        #receptor = self.get_coords(filename = receptor_filename)
        #rec_name = os.path.basename(receptor_filename)
        #rec_name = os.path.splitext(rec_name)[0]
        #rec_mol_name = os.path.splitext(rec_name)[0]
        # TODO this part will be handled by Fox

        # text data
        label = ["Absolute lowest energy", "Lowest energy in the largest cluster"]
        tag = ["LE", "LC"]
        time_info = time.localtime()[0:6]

        # pack the header:
        # start
        obuff  = "USER   ADVS_result> %d-%d-%d %d:%d:%d\n"  % (time_info)
        obuff += "USER   AD_rec> %s\n" % ( self.rec_mol_name)
        #obuff += "USER   AD_rec> %s\n" % ( receptor_filename)
        obuff += "USER   AD_runs,rmstol,tot_clusters> %d,%1.2f,%d\n" % ( total_runs, self.RMSTOL, len(clusters))
        obuff += "USER   AD_dlg_list> " # TODO think about a smart packing? # TODO convert to basename text

        # add the dlg list
        if not max_dlg_count or max_dlg_count > len(dlg_list):
            max_dlg_count = len(dlg_list)

        for d in range(max_dlg_count):
            obuff+="%s," % (os.path.basename( dlg_list[d]))
        obuff = obuff[:-1]+"\n"

        # add the results count  (1,2)
        obuff += "USER   AD_results> %d\n" % ( len(results) )

        # add the histogram
        obuff += "USER   AD_histogram> "
        for i in range(len(histogram)):
            obuff += "%2.3f:%d" % (histogram[i][0], histogram[i][1] )
            try:
                obuff += ":"+histogram[i][2]+","
            except:
                obuff += ","

        obuff = obuff[:-1]+"\n" # cut out the last ","
        #obuff += "REMARK AD_receptor> %s\n" % rec_name
        # separate pose MODELS are handled here
        for p in range(len(results)):
            # pose =  [L*_energy, L*_size, L*_pose] 
            pose = results[p]
            # find close contacts
            vdw_list = self.get_contact_atoms(pose[2], self.rec_crds)
            #vdw_list = self.get_contact_atoms(pose[2], receptor)
            # find hb 
            hb_dictionary = self.get_hb_interactions(pose[2]['text'], vdw_list)
            obuff += "MODEL    %d\n" % ( p+1)
            obuff += "USER   AD_label> %s\n" % (label[p])
            obuff += "USER   #     energy,\tleff,\tc_size,\tc_pc\n"
            obuff += "USER   AD_%s> %2.3f,\t%2.3f,\t%d,\t%3.2f\n" % (tag[p], results[p][0], self.lig_eff(pose[0], len(self.atype_list)), pose[1], self.percent(pose[1], total_runs) )

            # add hb acceptors info
            if hb_dictionary['acceptors']:
                obuff += "USER   AD_%s_hba> " % (tag[p])
                for pair in hb_dictionary['acceptors']:
                    obuff += "%s~~%s," % (self.pmv_atom_stripper(lig_mol_name, pair[0], contract = True), self.pmv_atom_stripper(self.rec_mol_name, pair[1],contract = True))
                obuff = obuff[:-1]+"\n"

            # add hb donors info
            if hb_dictionary['donors']:
                obuff += "USER   AD_%s_hbd> " % (tag[p])
                for pair in hb_dictionary['donors']:
                    obuff += "%s~~%s," % (self.pmv_atom_stripper(lig_mol_name, pair[0], contract = True), self.pmv_atom_stripper(self.rec_mol_name, pair[1], contract = True))
                obuff = obuff[:-1]+"\n"

            # add close contacts info
            if vdw_list:
                obuff += "USER   AD_%s_vdw> " % (tag[p])
                for a in vdw_list:
                    obuff += "%s," % (self.pmv_atom_stripper(self.rec_mol_name, a, contract = True))
                obuff = obuff[:-1]+"\n"
            for i in pose[2]['text']:
                if not i.startswith("MODEL") and not i.startswith("ENDMDL"):
                    obuff += i
            obuff +="ENDMDL   %d\n" % ( p+1 )
        return obuff



class DockingResultProcessor:
    """ Create pdbqt-plus file from a directory containing AutoDock results for a single ligand to a single receptor
        Input:
            directory containing one or more dlg files
        Options:
            rms_tolerance (default is 2.0)
            rms_refence_filename (default is to use input ligand coordinates)
            receptor_filename (default is set from gridmap names)
            best_energy_only (default is to create two pdbqt files:(1) best_energy (2) best_energy in largest_cluster)
            largest_cluster_only (default is to create two pdbqt files:(1) best_energy (2) best_energy in largest_cluster)
            largest_cluster_stem (default is 'ligandname_lc')
            max_cl_to_write (default is 10, use -1 to write them all)
            include_interactions (default is NOT to include interactions in pdbqt)
            detect_pi (default is NOT to detect pi-pi and cation-pi interactions)
            build_hbonds (default is to build hydrogen bonds)
            detect_close_contacts(default is to do so)

            
    """

    def __init__(self,  rms_tolerance=2.,
                        rms_reference=None,
                        receptor_filename=None,
                        write_both=True,
                        best_only=False,
                        largestCl_only=False,
                        le_stem="",
                        lc_stem="",
                        max_cl_to_write=-1,
                        include_interactions=False,
                        detect_pi=False,
                        build_hbonds=True,
                        detect_close_contacts=True,
                        intF=None,
                        docking=None
                        ):
        self.rms_tolerance=rms_tolerance
        self.rms_reference = rms_reference
        self.receptor_filename = receptor_filename
        self.write_both = write_both
        self.best_only = best_only
        self.largestCl_only = largestCl_only
        self.le_stem = le_stem
        self.lc_stem = lc_stem
        self.max_cl_to_write = max_cl_to_write
        self.include_interactions = include_interactions
        self.detect_pi = detect_pi
        self.build_hbonds = build_hbonds
        self.detect_close_contacts = detect_close_contacts
        self.intF = intF
        self.docking = docking
        if self.intF==None:
            self.intF = InteractionDetector(self.receptor_filename, detect_pi=self.detect_pi) #someday, detect_pi=True
        self.rmsTool = None
        if rms_reference is not None:
            self.rms_reference = rms_reference
            #file = os.path.join(directory, self.rms_reference)
            ref = Read(rms_reference)
            if not len(ref): 
                print "unable to read ", self.rms_reference
            else:
                ref = ref[0]
                coords = ref.allAtoms.coords
                refCoords = coords[:]
            self.rmsTool = RMSDCalculator(coords)

    
    def process(self, directory=None, outputfilename=None, verbose=False, rms_tolerance=None, dlg_list=[]):
        max_cl_to_write = self.max_cl_to_write
        #read all the docking logs in as one Docking
        if rms_tolerance is None:
            rms_tolerance=self.rms_tolerance
        if directory is not None:
            dlg_list = glob.glob(os.path.join(directory, "*.dlg"))
        if self.docking is None and  directory is None:
            print "either directory or dlg_list or docking must be specified"
            return "ERROR"
        if self.docking:
            d = self.docking
        else:
            d = Docking()
            ndlgs = 0
            for dlg in dlg_list:
                #code contributed by S.Forli
                try:  
                    d.readDlg(dlg)
                    ndlgs +=1
                except:
                    print "Error reading dlg file \"%s\""  % dlg
                    #break
        if not d.dlo_list:
            if self.docking:
                print "Sorry, provided docking not valid"
            else:
                print "Sorry, no valid DLG results found in the path : %s" % directory
            sys.exit(1)

        nconfs = len(d.ch.conformations)
        if verbose: print "processsing %d conformations from %d dlg files"%(nconfs, ndlgs)
        #if self.le_stem == "":
            #self.le_stem = origname+ "_vs_le"

        deduced_name = 0
        if self.receptor_filename is None:
            #must get it from the docking info
            for dpf_l in d.dlo_list[0].parser.dpfLines:
                if dpf_l.find('map')>-1:
                    break
            #'fld hsg1.maps.fld'.split() => ['fld', 'hsg1.maps.fld', '#', 'grid_data_file']
            try:
                rec_stem = dpf_l.split()[1].split('.')[0]
            except:
                print "exception deducing receptor_filename based on ", dpf_l
                exit(-1)
            self.receptor_filename = rec_stem + ".pdbqt"
            deduced_name = 1

        #check that it exists
        try:
            assert os.path.exists(self.receptor_filename)
        except:
            if deduced_name:
                print "Sorry unable to locate '%s' as deduced from gridmap names.\nPlease specify complete path to receptor file after '-r' "%( self.receptor_filename)
            else:
                print "Sorry unable to find receptor ", self.receptor_filename
            return

        rmsTool = self.rmsTool
        if rmsTool is None:
            #setup rmsd tool
            coords = d.ligMol.allAtoms.coords[:]
            rmsTool = RMSDCalculator(coords)
        d.clusterer.rmsTool = rmsTool
        d.clusterer.make_clustering(rms_tolerance) 
        #set up cl_dict here
        cl_dict = d.clusterer.clustering_dict[rms_tolerance]
        nclusters = len(cl_dict)
        cl_lengths = map(len, cl_dict)
        # adjust the number of clusters to write down if necessary
        if max_cl_to_write<0 or max_cl_to_write > nclusters:
            if verbose: print 'setting max_cl_to_write to nclusters ', nclusters
            max_cl_to_write = nclusters

        index_LC = cl_lengths.index(max(cl_lengths))
        largest_equals_best = True
        if index_LC!=0:
            largest_equals_best = False

        #ASSESS what to write:
        # both by default and if index_LC==0, ie largest is also best 
        # be only
        # lc only
        origname = d.ligMol.name
        if outputfilename is None:
            if largest_equals_best:
                d.ligMol.name = d.ligMol.name + "_vs"  #reserved for the best energy + largest cluster conf
            else:
                d.ligMol.name = d.ligMol.name + "_vs_le"  #reserved for the best energy 
            outputfilename = d.ligMol.name + ".pdbqt"

        ####-B best only
        ###best_only = False
        #be_outputfilename = d.ligMol.name + "_be.pdbqt"
        ####-L largestCl only
        ###largestCl_only = False
        # CORRECTION 6/29/2010
        # number of atoms used for ligand efficiency is NOW correctly heavy atoms ONLY!!
        total_num_lig_ats = len(d.ligMol.allAtoms)
        num_h_ats = len(d.ligMol.allAtoms.get(lambda x: x.element=='H'))
        num_lig_ats = total_num_lig_ats-num_h_ats
        #if self.write_both or self.best_only or index_LC==0:
        if self.write_both or self.best_only or index_LC==0:
            # create pdbqt file for the best
            best = cl_dict[0][0]
            # update the coordinates
            d.ch.set_conformation(best)
            new_coords = d.ligMol.allAtoms.coords
            #get index of this conformation
            best_ind = str(d.ch.conformations.index(best))
            # setup file for writing
            be_outfilename = os.path.join(directory , outputfilename)
            be_ptr = open(be_outfilename, 'w')
            #REMARK VirtualScreeningResult Date\n"
            sss, dr_str = d.clusterer.getInfoStr(comment='USER  AD> ',  ind=best_ind, rms=rms_tolerance,ncl_to_write=max_cl_to_write, include_dlgfilename_run=True)
            sss_list = sss.split('\n')
            dr_str_list = dr_str.split('\n') #dlgfilename and run info
            #look for omitted clusters
            omitted = ""
            if sss_list[0].find("omitted")>-1:
                omitted = sss_list[0] 
                sss_list = sss_list[1:]
            num_sss = len(sss_list)
            #print sss_list[0]
            if num_sss == 1:
                num_to_write = 1
            else:
                num_to_write = min(num_sss-1, max_cl_to_write)
            #for i in range(num_sss-1):  #sss starts with '#binding' so numbers are 1-based
            b_efficiency = best.binding_energy/num_lig_ats
            #be_ptr.write( "%s %4.2f %6.4f\n"%(sss_list[0],rms_tolerance, b_efficiency))
            #USER binding 2.0 200 runs 9 clusters
            be_ptr.write( "%s %4.2f %3d runs %3d clusters\n"%(sss_list[0],rms_tolerance, nconfs, nclusters))
            be_ptr.write( "USER  AD> ligand efficiency  %6.4f\n"%(b_efficiency))
            #be_ptr.write( "USER  AD> rmsd     LE    clu:size e_range  dlgfilename run ccrds? LE? LC?\n") #?CHANGE?
            be_ptr.write( "USER  AD> rmsd, LE, clu_size, clu_e_range, dlgfilename, run#, b_curCRDs, b_LE, b_LC\n") #?CHANGE?
            num_clusters = len(cl_dict)
            if num_clusters < num_to_write:
                num_to_write = num_clusters
            if verbose: print "num_to_write=", num_to_write
            for i in range(0, num_to_write):  #starts with '#binding' so numbers are 1-based
                cl = cl_dict[i]
                len_cl = len(cl)
                #add energy_range
                b_cl = cl[0]   #best in cluster
                w_cl = cl[-1]  #worst in cluster
                e_range = w_cl.binding_energy - b_cl.binding_energy
                if i==0: #BEST ENERGY
                    #eg:
                    #USER    0   0  0.000  0.000 -8.510  *
                    if index_LC!=0:
                        be_ptr.write( "%s,%d,%4.2f,%s,1,1,0\n"%(sss_list[i+1], len_cl,  e_range, dr_str_list[i]))
                    else:
                        #print "sss_list[", i,"] ", sss_list[i+1]
                        be_ptr.write( "%s,%d,%4.2f,%s,1,1,1\n"%(sss_list[i+1], len_cl,  e_range, dr_str_list[i]))
                elif i==index_LC: #Largest Cluster is NOT also lowest energy
                    be_ptr.write( "%s,%d,%4.2f,%s,0,0,1\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i]))
                #elif i==1: #correct first line USER   0.000 -5.800  *      1 0.00
                #    be_ptr.write( "%s,%d,%4.2f,%s,0,0,0\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i]))
                else:
                    be_ptr.write( "%s,%d,%4.2f,%s,0,0,0\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i]))
            if len(omitted): #if some clusters have been skipped
                be_ptr.write(omitted + "\n")
            #now write the rest of the ligand from its parser
            atm_ct = 0
            for line in d.ligMol.parser.allLines:
                #update the coordinates here...
                if line.find("HETATM")==0 or line.find("ATOM")==0:
                    nl = line[:30] +  "%8.3f%8.3f%8.3f"%(d.ligMol.allAtoms[atm_ct].coords[0], d.ligMol.allAtoms[atm_ct].coords[1], d.ligMol.allAtoms[atm_ct].coords[2])  + line[54:]
                    atm_ct += 1
                    #nl = line
                    be_ptr.write(nl + "\n")        
                else:    
                    be_ptr.write(line + "\n")
            be_ptr.close()
            if verbose: print " wrote best energy conformation to ", be_outfilename

            if index_LC==0 and verbose:
                print 'be == lc'
            # here's the point to include interaction information
            if verbose: print "include_interactions currently set to ", self.include_interactions
            if self.include_interactions:
                outputfilename = d.ligMol.name + ".pdbqt"
                #self.intF = InteractionDetector(self.receptor_filename, detect_pi=self.detect_pi)
                outputfilename = os.path.join(directory, outputfilename)
                if verbose: print "calling self.intF.processLigand with be_outfilename=", be_outfilename, "opf=", outputfilename
                sss2 = self.intF.processLigand(be_outfilename, outputfilename=outputfilename) 
            # add complete pathname here if not already in...
            if outputfilename==os.path.basename(outputfilename):
                outputfilename = os.path.join(directory, outputfilename)
            optr = open(outputfilename, 'r')
            lines = optr.readlines()
            optr.close()
            optr = open(outputfilename, 'w')
            optr.write("REMARK VirtualScreeningResult %s\n"%(time.asctime( time.localtime(time.time()))))
            for l in lines: 
                optr.write(l)
            optr.close()

        #   TODO #2: CLEAN THIS UP and test...
        if verbose:
            if index_LC==0: 
                print "LC is BE"
            print "index_LC =", index_LC
        if (self.write_both or self.largestCl_only) and index_LC!=0:
            # WRITE largestCL here!!
            #UPDATE  the ligand name
            # change the name so that correct strings are written in pdbqt   
            if self.lc_stem == "":
                self.lc_stem = origname+ "_vs_lc"
            d.ligMol.name = self.lc_stem
            outputfilename = d.ligMol.name + ".pdbqt"
            if verbose: 
                print "now d.ligMol.name = ", d.ligMol.name
                print "now outputfilename = ", outputfilename

            largest_cl = d.clusterer.clustering_dict[rms_tolerance][index_LC][0]
            d.ch.set_conformation(largest_cl)
            #sss += d.clusterer.getInfoStr(comment='USER ',  ind=index_LC, rms=rms_tolerance,ncl_to_write=max_cl_to_write)
            index_LCconf = d.ch.conformations.index(largest_cl)
            #sss, dr_str = d.clusterer.getInfoStr(comment='USER  AD> ',  ind=index_LCconf, rms=rms_tolerance,ncl_to_write=max_cl_to_write, include_dlgfilename_run=True)
            #check for omitted clusters
            #omitted = ""
            #sss_list = sss.split('\n')
            #dr_str_list = dr_str.split('\n')
            #if sss_list[0].find("omitted")>-1:
            #    omitted = sss_list[0]
            #    sss_list = sss_list[1:]
            #num_sss = len(sss_list)
            
            if num_sss == 1:
                num_to_write = 1
            else:
                num_to_write = min(num_sss-1, max_cl_to_write)
            #sss_list = sss.split('\n')
            # open the lc file
            lc_outfilename = os.path.join(directory, self.lc_stem +".pdbqt")
            lc_ptr = open(lc_outfilename, 'w')
            #lc_ptr.write( "%s %4.2f\n"%(sss_list[0],rms_tolerance))
            #USER binding 2.0 200 runs 17 clusters 
            lc_ptr.write( "%s %4.2f %3d runs %3d clusters\n"%(sss_list[0],rms_tolerance, nconfs, nclusters))
            lc_b_efficiency = largest_cl.binding_energy/num_lig_ats
            lc_ptr.write( "USER  AD> ligand efficiency  %6.4f\n"%(lc_b_efficiency))
            #lc_ptr.write( "USER  AD> rmsd     LE_LC  cl:size e_range\n")
            #lc_ptr.write( "USER  AD> rmsd     LE_LC    clu:size e_range  dlgfilename run\n")
            lc_ptr.write( "USER  AD> rmsd, LE_LC, clu_size, clu_e_range, dlgfilename, run#, b_curCRDs, b_LE, b_LC\n") #?CHANGE?
            for i in range(0, num_to_write):  #sss_list starts with '#binding' so indicies into it are 1-based
                cl = d.clusterer.clustering_dict[rms_tolerance][i]
                len_cl = len(cl)
                #add energy_range
                b_cl = cl[0]
                w_cl = cl[-1]
                e_range = w_cl.binding_energy - b_cl.binding_energy
                if i==0: # index_LC==0 handled with BE above
                    lc_ptr.write( "%s,%d,%4.2f,%s,0,1,0\n"%(sss_list[i+1], len_cl,  e_range, dr_str_list[i]))
                elif i==index_LC:
                    m_cl = d.clusterer.clustering_dict[rms_tolerance][index_LC]
                    mb_cl = m_cl[0]
                    mw_cl = m_cl[-1]
                    len_cl = len(m_cl)
                    lc_ptr.write( "%s,%d,%4.2f,%s,1,0,1\n"%(sss_list[i+1][:-1], len_cl, e_range, dr_str_list[i]))
                else:
                    lc_ptr.write( "%s,%d,%4.2f,%s,0,0,0\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i]))
            if len(omitted): #if some clusters have been skipped
                lc_ptr.write(omitted + "\n")
            #now write the rest of the ligand from its parser
            atm_ct = 0
            for line in d.ligMol.parser.allLines:
                if line.find("HETATM")==0 or line.find("ATOM")==0:
                    nl = line[:30] +  "%8.3f%8.3f%8.3f"%(d.ligMol.allAtoms[atm_ct].coords[0], d.ligMol.allAtoms[atm_ct].coords[1], d.ligMol.allAtoms[atm_ct].coords[2])  + line[54:]
                    lc_ptr.write(nl + "\n")        
                    atm_ct += 1
                else:
                    lc_ptr.write(line+"\n")
            lc_ptr.close()
            #print " wrote lc to ", lc_outfilename
            if verbose: print " wrote lc to ", lc_outfilename
            if self.include_interactions:
                #intF = InteractionDetector(self.receptor_filename, self.detect_pi) #someday, detect_pi=True
                outputfilename = os.path.join(directory, outputfilename)
                if verbose: print "calling self.intF.processLigand with lc_outfilename=", lc_outfilename, "opf=", outputfilename
                sss2 = self.intF.processLigand(lc_outfilename, outputfilename=outputfilename) 
            if outputfilename==os.path.basename(outputfilename):
                outputfilename = os.path.join(directory, outputfilename)
            optr = open(outputfilename, 'r')
            lines = optr.readlines()
            optr.close()
            optr = open(outputfilename, 'w')
            optr.write("REMARK VirtualScreeningResult %s\n"%(time.asctime( time.localtime(time.time()))))
            for l in lines: optr.write(l)
            optr.close()


class VinaResultProcessor:
    """ Creates pdbqt-plus files from AutoDock Vina docking results to a single receptor
        Methods:
            __init__
                Required arguments:
                    receptor_filename,
                Optional arguments:
                    best_only -only output best conformation(default is to create a pdbqt files for each MODEL)
                    max_num_to_write (default is to write one pdbqt file for each MODEL)
                    detect_close_contacts (default is to do so)
                    detect_pi (default is NOT to detect pi-pi and cation-pi interactions)
                    intV (default is to use default InteractionDetector)
            
            process 
                Required argument:
                    AutoDockVina resultfile containing 1 or more MODELS
                Optional arguments:
                    outputfilestem (default is '_Vvs')
                    output_stem (default is to use ligname + '_vvs'+index)
                    verbose
    """

    def __init__(self,  receptor_filename,
                        best_only=False,
                        max_num_to_write=-1,
                        detect_close_contacts=True,
                        detect_pi=False,
                        intD=None,
                        ):
        #@@ receptor_filename is required !!
        self.receptor_filename = receptor_filename
        #check that receptor_filename exists
        try:
            assert os.path.exists(receptor_filename)
        except:
            print "Sorry unable to find receptor file ", receptor_filename
            exit()
        self.best_only = best_only
        self.max_num_to_write = max_num_to_write
        if self.best_only: 
            self.max_num_to_write = 1
        self.detect_pi = detect_pi
        self.detect_close_contacts = detect_close_contacts
        self.intD = intD
        if self.intD==None:
            self.intD = InteractionDetector(self.receptor_filename, detect_pi=self.detect_pi) 

    
    def process(self, vina_resultfile, outputfilestem=None, remove_model_str=False, best_only=None, verbose=False):
        if best_only is None:
            best_only = self.best_only
        ligMols = Read(vina_resultfile) 
        #NB:
        #by default Read returns multiple molecules with different sets of coords
        #modelsAs='conformations' returns 1 molecule with multiple sets of coords
        if not len(ligMols):
            print "unable to read vina resultfile ", vina_resultfile
            sys.exit()
        #nconfs
        ligname = ligMols[0].name
        nconfs = len(ligMols)
        if verbose: print "%s contains %d conformations for ligand %s"%(vina_resultfile, nconfs, ligname)
        num_to_write = self.max_num_to_write
        if best_only:
            num_to_write=1
        elif num_to_write < 0 or nconfs < num_to_write:
            num_to_write = nconfs
        #set up best
        best = ligMols[0]
        parser_lines = best.parser.allLines
        receptor_filename = self.receptor_filename
        ligname, ext_type = os.path.splitext(os.path.basename(best.parser.filename))
        vina_results = best.vina_results
        # eg: list of vina results attached to first MODEL 'best'
        #[['-11.6', '0.000', '0.000'], ['-10.6', '2.048', '11.070'], 
        # ['-10.6', '1.360', '4.236'], ['-10.2', '1.871', '10.935'], 
        # ['-10.1', '2.121', '10.498'], ['-10.0', '1.473', '2.408'], 
        # ['-9.9', '1.856', '10.886'], ['-9.9', '2.506', '11.224'], 
        # ['-9.9', '2.366', '10.659']]
        # number of atoms used for ligand efficiency is NOW correctly heavy atoms ONLY!!
        total_num_lig_ats = len(best.allAtoms)
        num_h_ats = len(best.allAtoms.get(lambda x: x.element=='H'))
        num_lig_ats = total_num_lig_ats - num_h_ats
        if verbose: print "processsing %d conformations from vina result %s file"%(nconfs, vina_resultfile)
        start_str = "USER  AD> "
        vsV_str = "REMARK AutoDockVina VirtualScreeningResult "
        if verbose: print "num_to_write=", num_to_write
        #11/2010:previously next line was "for i in range(len(ligMols)):"
        for i in range(num_to_write):
            curLig = ligMols[i]
            curName = curLig.name
            if remove_model_str:
                modelIndex=curName.find("_model")
                if modelIndex>-1:
                    curName = curName[:modelIndex]+ curName[modelIndex+6:]
                    if verbose: print "now curName=", curName
                curLig.name = curName
            #setup outputfilestem
            if outputfilestem is None:
                outputfilename = curName +".pdbqt"
            else:
                outputfilename = outputfilestem+"%s"%str(i+1) + ".pdbqt"
                #outputfilename = curName +outputfilestem+"%s"%str(i+1) + ".pdbqt"
            # create pdbqt+ file for this docked conf
            #1. open file for writing
            #outfilename = os.path.join(directory , outputfilename)
            optr = open(outputfilename, 'w')
            #2. write the ligand using its parser's lines
            atm_ct = 0
            found_model = 0
            for ctr in range(len(parser_lines)):
                line = parser_lines[ctr]
                if line.find("MODEL")==0:
                    found_model += 1
                    if found_model>1:
                        if verbose: print "found ", i,"-TH MODEL!",
                        parser_lines = parser_lines[ctr:]
                        if verbose: print "NOW len(parser_lines)=", len(parser_lines), '\n'
                        break
                optr.write(line)
            if verbose: print " wrote %dth energy conformation to %s"%(i, outputfilename)
            optr.close()
            remove_modelStr = found_model and remove_model_str
            self.intD.processLigand(outputfilename, outputfilename=outputfilename, buildHB=0, remove_modelStr=remove_modelStr)
            optr = open(outputfilename, 'r')
            lines = optr.readlines()
            optr.close()
            #4. write header HERE instead
            optr = open(outputfilename, 'w')
            ostr = "%s %s\n"%(vsV_str, time.ctime())
            optr.write(ostr)
            optr.write( "%s %d of %d MODELS\n"%(start_str, i+1, nconfs ))
            b_efficiency = curLig.vina_energy/num_lig_ats
            optr.write( "%s ligand efficiency  %6.4f\n"%(start_str, b_efficiency))
            energy, min_rms, max_rms = map(float, best.vina_results[i])
            optr.write( "%s %4.2f, %6.3f, %6.3f\n"%(start_str, energy, min_rms, max_rms))
            for l in lines:
                optr.write(l)
            optr.close()
            


