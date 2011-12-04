## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2002
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/DlgParser.py,v 1.60.2.4 2009/09/17 19:54:59 rhuey Exp $
#
# $Id: DlgParser.py,v 1.60.2.4 2009/09/17 19:54:59 rhuey Exp $
#
#
#
#
#
#
#

"""
This Object parses the result of an AutoDock job and returns a dictionary. 

"""
import os
from string import find, join, replace, split, rfind, strip
import re
import numpy.oldnumeric as Numeric

from AutoDockTools.ResultParser import ResultParser


class DlgParser(ResultParser):
    """ reads log from a AutoDock docking and return structured data"""

    keywords = ResultParser.keywords + [
        'coords',
        'vdw_energies',
        'estat_energies',
        'total_energies',   #vdw_energies+estat_energies
        'inhib_constant',
        'intermol_energy',  #(1) 
        'internal_energy',  #(2) NB: 1+2->final docked energy
        'torsional_energy', #(3) NB: 1+3->free energy of binding
        'run',
        'parameter_file',   #AD4 parameter library filename
        'include_1_4_interactions',   #AD4 internal energy switch
        ]
        

    def __init__(self, dlgFile=None):
        """selected dlgFile,ok sets which docked conformations to show"""
        ResultParser.__init__(self)
        self.filename = dlgFile
        self.WARNINGS = []
        self.population_list = []
        self.population = []
        #set up dict here
        self.getReDict()
        self.ligand_atom_count = 0
        if dlgFile:
            self.filename = os.path.basename(dlgFile)
            self.parse(dlgFile)



    def parse(self, filename):
        """
        uses keys:
            '|           AutoDock' 
            'DPF> outlev'
            'DPF> '
            'INPUT-PDBQ: '
            'DOCKED: '
            'State='
            '^Seeds'
            '^   [1-9]
        after parsing: 
            self.outlev is switch to correctly find coord fields
            self.clist is list of dictionaries for states,
            self.dpfLines are for building DockingParameters object
            self.histogramlines show results histogram
            self.clusterlines have clustering info
            self.modelList is list of dicts of docked: coords plus energy info 
            self.wroteAll is a flag set by 'write_all_cluster_members' or
                    by 'write_all'
        """
        self.filename = filename
        self.version = 4.2
        #reset
        for item in ['clist','clusterlines','dpfLines','histogramlines',\
                     'ligLines','modelList', 'initial_clist']:
            setattr(self, item, [])

        self.wroteAll = 0
        self.clusterRecord = None
        self.getReDict()
        try:
            dlgptr = open(filename, 'r')
        except:
            raise IOError
        allLines = dlgptr.readlines()
        self.allLines = allLines
        self.match(allLines)
        if self.version>=4.0:
            ad4_keywds = ['vdw_energy', 'estat_energy', 
                          'inhib_constant_units',
                          'vdw_hb_desolv_energy',
                          'electrostatic_energy',
                          'moving_ligand_fixed_receptor',
                          'moving_ligand_moving_receptor',
                          'total_internal',
                          'ligand_internal',
                          'receptor_internal',
                          'unbound_energy',
                          ]
            self.keywords.extend(ad4_keywds)


    def getReDict(self):
        if hasattr(self, 'reDict'):
            for k, d in self.reDict.items():
                d['lines'] = []
            return
        self.reDict = {}
        #had to add MODEL|^USER|^ATOM key for outlev -1
        #seems to slow stuff down alot
        self.compute_unbound_extended = False
        self.reKeys = [
            '                 |            AutoDock', #used to correct autodock4 output 
            'DPF> outlev',       #format switch
            'DPF> compute_unbound_extended', # affects number of seeds
            'DPF> write_all',    #wroteAll switch
            'INPUT-PDBQ: USER    NEWDPF',         #extra info in case of clustering dlg
            'DPF> ',             #lines for dpo
            'INPUT-PDBQ: ',      #lines for ligand
            'INPUT-PDBQT: ',     #lines for ligand, v4.0_old_version
            'INPUT-LIGAND-PDBQT: ',     #lines for ligand, v4.0
            'INPUT-FLEXRES-PDBQT: ',     #lines for flex_res, v4.0
            'DOCKED: ',          #lines for models
            '^MODEL|^USER|^ATOM|^ENDMDL',       #outlev -1 modellines
            'State=',            #lines for conformations
            '^ [1-9]|^  [1-9]|^   [1-9]',   #cluster and histogram lines
            '^Seeds',            #lines for seed
            'Total number of atoms found in PDBQT file', #to get ligand at #
            '^Atom: ID',         #for non-bond table
            'WARNING',           #accumulate all warning messages
            'Energy= ',          #accumulate initial states available if outlev>1
            '<population size=',
            ".*autodock.*: Successful Completion",
            "Coordinates of Central Grid Point of Maps",
            ]
        self.reFuncs = [ 
            self.setADVersion,
            self.setOutlev,
            self.setComputeUnboundExtended,
            self.setWroteAll,
            self.getNewDpfInfo,
            self.processDpfLines,
            self.processLigLines,
            self.processLigLinesV4,
            self.processLigLinesV4,
            self.processFlexResLinesV4,
            self.getModelLines,
            self.getShortModelLines,
            self.getDlgStates,
            self.getClusterInfo,
            self.getSeedInfo,
            self.getLigandAtomCount,
            self.getNonBondTable,
            self.getWARNINGS,
            self.getInitialConfs,
            self.getPopulations,
            self.getTime,
            self.getCenterMapPt,
            ]

        for i in range(len(self.reKeys)):
            k = self.reKeys[i]
            dict =  self.reDict[k] = {}
            dict['re'] = re.compile(k)
            dict['lines'] = []
            dict['func'] = self.reFuncs[i]


    def match(self, allLines, verbose=False):
        self.getReDict()
        self.tested = 1
        for i in range(len(allLines)):
            item = allLines[i]
            #if find item, mark it found + don't test
            for k in self.reKeys:
                d = self.reDict[k]
                m = d['re'].match(item)
                if m:
                    d['lines'].append(item) 
                    break
            #process Energy= separately
            #energy_key = 'Energy= '
            energy_key = "Energy= "
            d = self.reDict[energy_key]
            if item.find(energy_key)>-1:
                d['lines'].append(item)
                d['lines'].append(allLines[i+1])
        for k in self.reKeys:
            d = self.reDict[k]
            lines = d['lines']
            if lines==[]:
                input_key = 'INPUT-PDBQ: '
                if k==input_key:
                    if verbose and self.version<4.0:
                        print "!!no lines found for key=", k, "!!"
                else:
                    if verbose:
                        print "!!no lines found for key=", k, "!!"
            else:
                apply(d['func'], (lines,), {})


    def setADVersion(self,lines):
        if len(lines):
            for l in lines:
                if find(l, 'AutoDock')>-1:
                    ll = split(l)
                    version = ll[2]
                    if len(version)>3:
                        version = version[:3]
                    self.version = float(version)
                    break
        else:
            self.version = 3.0


    def setOutlev(self,lines):
        if len(lines):
            ll = split(lines[0])
            self.outlev = int(ll[2])
        else:
            self.outlev = None


    def setComputeUnboundExtended(self,lines):
        if len(lines):
            self.compute_unbound_extended = True


    def setWroteAll(self, lines):
        if len(lines):
            self.wroteAll = 1


    def getLigandAtomCount(self, lines):
        self.ligand_atom_count = int(split(lines[0])[-2])
        #print "parsed ligand_atom_count =", self.ligand_atom_count


    def getWARNINGS(self, lines):
        #print "found ", len(lines), ' WARNINGS'
        self.WARNINGS = lines


    def getInitialConfs(self, lines):
        #print "found ", len(lines), ' InitialConfs'
        self.InitialConfs= lines
        x = len(lines)/2
        for j in range(x):
            energy = float(lines[2*j].strip().split()[-1])
            clist = map(float, lines[2*j+1].strip().split())
            self.initial_clist.append([energy, clist])


    def getCenterMapPt(self, lines):
        for l in lines:
            if l.find('Coordinates of Central Grid Point of Maps')>-1:
                ll = split(l, '=')
                str_pt = ll[1].strip()  #"(13.231, 56.852, 54.080)
                #str_pt = ll[1]  #"(13.231, 56.852, 54.080)
                self.center_pt = map(float, str_pt[2:-2].split(','))


    def getTime(self, lines, echo=False):
        ind = self.allLines.index(lines[-1])
        #print "ind=", ind
        self.info_line = None
        self.total_time = 0
        for info_line in self.allLines[ind+1:]:
            #print "testing info_line->", info_line
            #print "info_line.find('Real=')=", info_line.find('Real=')
            if info_line.find("Real=")>-1:
                #print "FOUND IT!!"
                self.info_line = info_line
                break
        if self.info_line:
            #Real= 23m 38.38s,  CPU= 22m 08.78s,  System= 0.49s
            time_string = self.info_line.split(',')[0][6:]
            #23m 38.38s
            time_list = time_string.split()
            #['5h', '14m', '39.77s']
            #['14m', '39.77s']
            #['39.77s']
            scale_factor ={ 'h':60*60,'m':60,'s':1}
            for x in time_list:
                #['5h', '14m', '39.77s']
                self.total_time += scale_factor[x[-1]] * float(x[:-1])
        

    def getPopulations(self, lines, echo=False):
        currently_build_pop = True
        ind = self.allLines.index(lines[0])
        if echo: print "found start of population at index=", ind
        found = False
        all_population_lines = []
        plist = []
        currently_building_pop = True
        for line in self.allLines[ind:]:
            if line.find("</pop")==0:
                if echo: print 'found end of a population', ind
                found = True
                all_population_lines.append(plist)
                if echo: print "added plist:", len(plist)
                currently_building_pop = False
                plist = []
                #break
            elif line.find("<pop")==0:
                if echo: print 'found start of a population', ind
                currently_building_pop = True
                plist = []
            elif currently_building_pop:
                plist.append(line)
        x = len(plist)
        if echo: print x, 'lines in getInitialPopulation'
        self.population_list = []
        #figure out number of torsions:
        num_values = len(all_population_lines[0][0].strip().split())
        # 0       1       2   3      4      5      6      7     8     9 # ....
        #index   enrgy   age  trnX  trnY   trnZ    qtX    qtY   qtZ  qtANG [tor1 tor2 ...]
        #150	3.36e+05  0	 40.996 13.492 23.851  0.092 -0.658 0.748 9.890  
        num_torsions = num_values - 10
        for plist in all_population_lines:
            #previously
            #150	3.36e+05	40.996 13.492 23.851  0.092 -0.658 0.748 9.890  
            #current ad4 version has age between energy and tran_x, eg 0, third number
            #150	3.36e+05  0	 40.996 13.492 23.851  0.092 -0.658 0.748 9.890  
            current_pop = []
            for line in plist:
                ll = line.strip().split()
                ind_index = int(ll[0])
                age  = int(ll[2])
                clist = map(float, ll[3:])
                cdict = {}
                cdict['binding_energy'] = float(ll[1])
                cdict['age'] = int(ll[2])
                cdict['trn_x'] = float(ll[3])
                cdict['trn_y'] = float(ll[4])
                cdict['trn_z'] = float(ll[5])
                cdict['qtn_nx'] = float(ll[6])
                cdict['qtn_ny'] = float(ll[7])
                cdict['qtn_nz'] = float(ll[8])
                cdict['qtn_ang_deg'] = float(ll[9])
                cdict['num_torsions'] = num_torsions
                if num_torsions>0:
                    cdict['torsion_values'] = clist[7:]
                current_pop.append(cdict)
            self.population_list.append(current_pop)
        self.initial_pop = self.population_list[0]
        self.population = self.initial_pop
            

    def getNonBondTable(self, lines, echo=False):
        ind = self.allLines.index(lines[0])
        ct = self.ligand_atom_count
        nb_lines = self.allLines[ind+2:ind+2+ct]
        self.nb_array = Numeric.zeros((ct, ct))
        if echo:
            for l in nb_lines:
                print l,
        for i in range(ct):
            l = nb_lines[i][10:-1]
            for j in range(ct):
                jind = j*2 + 1
                if l[jind]=='X':
                    self.nb_array[i][j]=1
            if echo:
                print
                print l
                print self.nb_array[i]


    def getSeedInfo(self, lines):
        seeds = []
        for l in lines:
            ilist = split(l)
            slist = []
            for item in ilist[1:]:
                slist.append(int(item))
            seeds.append(tuple(slist))
            #seeds.append(tuple((int(ilist[1]),int(ilist[2]))))
        if self.version>=4.0 and self.compute_unbound_extended:
            #there are seeds for compute_unbound_extended
            #remove the first entry
            seeds = seeds[1:]
        if len(seeds)==len(self.clist):
            for i in range(len(seeds)):
                conf = self.clist[i]
                if conf.has_key('run') or self.version>=4.0:
                    #outlev -1 prints States in order
                    s = seeds[i]
                    conf['rseed1'] = s[0]
                    if len(s)==2:
                        conf['rseed2'] = s[1]
                    #elif not conf.has_key('run'):
                #else:
                #    continue
                    #outlev -1 prints States in order
                    #s = seeds[i]
                    #conf['rseed1'] = s[0]
                    #if len(s)==2:
                        #conf['rseed2'] = s[1]
                #else:
                    #run = conf['run']
                    ##run is 1-based; seeds indexed starting w/0
                    #s = seeds[run-1]
                    #conf['rseed1'] = s[0]
                    #if len(s)==2:
                        #conf['rseed2'] = s[1]
                #else:
                #run = conf['run']
                ##run is 1-based; seeds indexed starting w/0
                #s = seeds[run-1]
                #conf['rseed1'] = s[0]
                #if len(s)==2:
                #    conf['rseed2'] = s[1]


    def getClusterInfo(self, lines):
        cl = self.clusterlines
        hl = self.histogramlines
        for l in lines:
            if find(l, 'RANKING')>-1:
                cl.append(l[:-1])
            elif find(l, '#')>-1:
                hl.append(l[:-1])
        if len(cl):
            self.getClusterRecord(cl)
        else:
            self.clusterRecord = None


    def getClusterRecord(self, cl):
        #print 'in getClusterRecord'
        clRecList = []
        curList = []
        #curList gets list of conf info
        curInd = int(split(cl[0])[0])
        ctr = 1
        for l in cl:
            ll = split(l)
            #when built, newList is
            #[Rank,SubRank,Run,DockedEnergy,ClusterRMSD,RefREMSD]
            newList = map(lambda x:int(x),ll[:3])
            #3/29/05
            if self.wroteAll and self.version<4.0:
                #print "setting run number to ", ctr
                newList[2] = ctr
                ctr = ctr + 1
            newList2 = map(lambda x:float(x),ll[3:-1])
            newList.extend(newList2)
            if newList[0]==curInd:
                curList.append(newList)
            else:
                clRecList.append(curList)
                curList = [newList]
                curInd = newList[0]
        clRecList.append(curList)
        self.clusterRecord = clRecList
        
            

    def processDpfLines(self, lines):
        for l in lines:
            if len(l[:-1])>5:
                #check for any garbage in dpf
                if l.find('DPF> #')==0:
                    continue 
                if l.find('move')>-1:
                    #check for run_on ligand_filename here
                    if l.find('pdbqt#')>-1:
                        ind = l.find('pdbqt#')
                        l = l[:ind+5] + ' ' + l[ind+5:]
                self.dpfLines.append(l[5:-1])
                if l.find('ga_run')>-1:
                    ll = split(l)
                    self.runs = int(ll[2])
                    #print "self.runs=", self.runs
                elif l.find('runs')>-1 and l.find('simanneal')<0:
                    ll = split(l)
                    self.runs = int(ll[2])
                    #print "SA: self.runs=", self.runs

    
    def getNewDpfInfo(self, lines):
        if not len(lines):
            return
        dpfLines = []
        keys = []
        for l in lines:
            ll = split(l)
            ind = ll.index('NEWDPF')
            k = ll[ind+1]
            #only add each key once
            if k not in keys:
                dpfLines.append(join(ll[ind+1:]))
                keys.append(k)
        if len(dpfLines):
            self.dpfLines.extend(dpfLines)


    def processLigLines(self, lines):
        #do not use this for version 4.0
        if self.version>=4.0:
            return
        ligLINES = []
        foundRun = 0
        for l in lines:
            #in clustering dlg, multiple copies of input-pdbq are present
            if find(l, ' Run')>-1 and foundRun:
                break
            elif find(l, ' Run')>-1:
                foundRun = 1
            else:
                ligLINES.append(l[12:-1])
        #check here to remove lines of just spaces
        nl = []
        for l in ligLINES:
            if len(strip(l)):
                nl.append(l)
        self.ligLines = nl
        #print "in processLigLines:len(ligLines)=", len(nl)
        #self.ligLines = ligLINES


    def processFlexResLinesV4(self, lines):
        #print "in processFlexResLinesV4: len(self.ligLines=)", len(self.ligLines)
        if self.version<4.0:
            print "not version 4.0! RETURNING!!"
            return
        ligLINES = []
        foundRun = 0
        ind = 21
        for l in lines:
            #in clustering dlg, multiple copies of input-pdbq are present
            if find(l, 'Run')>-1 and foundRun:
                break
            elif find(l, 'Run')>-1:
                foundRun = 1
            elif find(l, '^_____________________')>-1:
                #last line is ________________-
                break
            else:
                ligLINES.append(l[ind:-1])
        #check here to remove lines of just spaces
        nl = []
        for l in ligLINES:
            if len(strip(l)):
                nl.append(l)
        self.flex_res_lines = nl
        #print "end pFRLV4: len(self.flex_res_lines)=", len(nl)
        #print "end processFlexResLinesV4: len(self.ligLines=)", len(self.ligLines)
        self.hasFlexRes = True
        self.flex_res_count = nl.count("REMARK  status: ('A' for Active; 'I' for Inactive)")

    
    def processLigLinesV4(self, lines):
        #use this only for versions 4.0 and greater
        #hack to fix run-on INPUT-PDBQT: line from wcg 
        first_line = lines[0]
        if len(first_line)>100 and first_line.count("INPUT-PDBQT:")>2:
            lines = first_line.replace("INPUT-PDBQT:", "\nINPUT-PDBQT:").split('\n')
            #exclude leading '' and trailing '' from split
            lines = lines[1:-1]
            #print "replaced lines with ", len(lines), ' lines'
        if self.version<4.0:
            #print "not version 4.0"
            return
        ligLINES = []
        foundRun = 0
        if find(lines[0], 'INPUT-LIGAND-PDBQT')==0:
            ind = 20
        elif find(lines[0], 'INPUT-PDBQ')==0:
            ind = 13

        for l in lines:
            #in clustering dlg, multiple copies of input-pdbq are present
            if find(l, 'Run')>-1 and foundRun:
                break
            elif find(l, 'Run')>-1:
                foundRun = 1
            elif find(l, 'TORSDOF')>-1:
                #eg run-on TORSDOF line: 
                #TORSDOF 3___
                lastChar = l.find('_')
                #previously l[13:-1]
                #ligLINES.append(l[13:lastChar])
                ligLINES.append(l[ind:lastChar])
                l_index = lines.index(l)
                if l_index==len(lines)-1:
                    #print "found TORSDOF on last line!"
                    break
                next_line = lines[l_index+1]
                if find(next_line, 'BEGIN_RES')<0:
                    print len(lines[l_index:]), ' lines left unparsed!!!'
                    break
            else:
                #ligLINES.append(l[13:-1])
                ligLINES.append(l[ind:-1])
        #check here to remove lines of just spaces
        nl = []
        for l in ligLINES:
            if len(strip(l)):
                nl.append(l)
        #print "end pLV4: len(self.ligLines)=", len(nl)
        self.ligLines = nl


    def getModelLines(self, lines):
        #print "in getModelLines with ", len(lines), ' lines'
        if not len(lines):
            return
        modelList = []
        if find(lines[0], 'DOCKED')==0:
            ind = 8
        elif find(lines[0], 'INPUT-LIGAND-PDBQT')==0:
            ind = 20
        elif find(lines[0], 'INPUT-PDBQ')==0:
            ind = 12
        #preprocess lines to remove DOCKED or INPUT-PDBQ
        nlines = []
        ctr = 0
        self.run_models = {}
        has_docked = False
        for l in lines:
            #if self.version!=4.0:
            #    nlines.append(l[ind:-1])
            #else:
            #    if find(l, 'ATOM')>-1:
            #        newLine = l[ind:64] + l[66:-1]
            #        nlines.append(newLine)
            #    else:
            #        nlines.append(l[ind:-1])
            if self.version>=4.0:
                if l.find("DOCKED: USER    Run = ")>-1:
                    ll = l.split()
                    ctr = int(ll[4])
                    #print "found run ", ctr
                    has_docked = True
                    nlines = []
                elif l.find("DOCKED: ENDMDL")>-1:
                    self.run_models[ctr] = nlines
                    #print "saved run ", ctr, " nlines=", len(nlines)
                    has_docked = False
                    ctr = -1
                elif has_docked is True:
                    #print 'appending ', l
                    if l.find("ATOM")>-1:
                        try:
                            nlines.append(l[ind:-1])
                        except:
                            print 'cutting out 64+65 of ', l
                            nlines.append(l[ind:64]+ l[66:-1])
                    else:
                        nlines.append(l[ind:-1])
            else:
                nlines.append(l[ind:-1])
        #print "len(nlines)=", len(nlines)
        if self.version<4.0:
            curMod = [nlines[0]]
            for l in nlines[1:]:
                if find(l, 'MODEL')>-1:
                    modelList.append(curMod)
                    curMod = [l]
                else:
                    curMod.append(l)
            modelList.append(curMod)
            #self.makeModels(modelList)
        else:
            #print "setting modelList to ", self.run_models.values()
            modelList = []
            for i in range(self.runs):
                try:
                    modelList.append(self.run_models[i+1])
                except:
                    print "unable to load expected model number ", i+1
        self.modelList = modelList
        self.makeModels(modelList)


    def getShortModelLines(self, lines):
        if not len(lines):
            return
        modelList = []
        curMod = [lines[0][:-1]]
        endmdl = 0
        for l in lines[1:]:
            #use endmdl because cluster dlg format
            # has MODEL then lines w/o MODEL then MODEL
            #again in desc of 1 model...
            if find(l, 'ENDMDL')>-1:
                curMod.append(l)
                modelList.append(curMod)
                endmdl = 1
            elif endmdl:
                curMod = [l[:-1]]
                endmdl = 0
            else:
                curMod.append(l[:-1])
        #don't append after for loop 
        #because it happens in the else
        #print 'len(modelList)=', len(modelList)
        self.modelList = modelList
        self.makeModels(modelList)


    def getDlgStates(self, lines):
        #print "in getDlgStates: len(self.clist)=", len(self.clist)
        if len(self.reDict['^MODEL|^USER|^ATOM|^ENDMDL']['lines']):
            #print 'not building states because models present'
            return
        if len(self.clist)==self.runs:
            for i in range(len(self.clist)):
                cl = self.clist[i]
                if not cl.has_key('run'):
                    cl['run'] = i+1
            return
        else:
            #print "resetting self.clist
            self.clist = []
        for l in lines:
            # in test-1 State= + 17 items: 3 trans, 4quat + ndihe(10) torsions
            xx = split(l)
            # remove possible punctuation
            for ind in range(len(xx)):
                if xx[ind][-1]==',': 
                    xx[ind] = xx[ind][:-1]
                    raise 'comma'
                if xx[ind][-1]=='.': 
                    xx[ind] = xx[ind][:-1]
                    raise 'period'
                
            #transList = xx[1:4]
            trans = []
            for p in [0,1,2]:
                trans.append(float(xx[p+1]))
            trans = tuple(trans)

            quat = []
            for n in xx[4:8]:
                quat.append(float(n))
            quat = tuple(quat)

            angList = []
            #NB: here torsions are in the same line
            for n in xx[8:]:
                angList.append(float(n))

            #BUILD A DICTIONARY and put it in clist
            d = {}
            d['trn_x'] = trans[0]
            d['trn_y'] = trans[1]
            d['trn_z'] = trans[2]
            d['qtn_nx'] = quat[0]
            d['qtn_ny'] = quat[1]
            d['qtn_nz'] = quat[2]
            d['qtn_ang_deg'] = quat[3]
            d['num_torsions'] = len(angList)
            d['torsion_values'] = angList
            self.clist.append(d)
            #print "added ", len(self.clist), "th conformation"


    def makeModels(self, modelList):
        #print "in makeModels with ", len(modelList), ' models to build'
        #print "self.version=", self.version
        if hasattr(self, 'clist') and len(self.clist):
            #print "parser already has clist, returning!"
            return
        #print "##########  SETTING clist to [] ###########"
        clist = []
        ctr = 1
        for curMod in modelList:
            clist.append(self.makeModel(curMod))
            #@@@CHECK THIS@@@
            clist[-1]['run'] = ctr
            ctr +=1
            #print "now len(clist)=", len(clist)
        self.clist = clist
        #print "end of  makeModels: len(self.clist)=", len(self.clist)


    def makeModel(self, lines):
        coords = []
        vdW = []
        Elec = []
        d = {}
        corr = 0
        if self.outlev==-1:
            corr = -1
        binding_energy2 = None
        version = self.version
        for l in lines:
            ll = split(l)
            #if find(l, 'MODEL')>-1:
            #    d['num'] = int((ll)[1])
            if find(l, 'Run')>-1 and find(l, 'Rank')==-1:
                d['run'] = int((ll)[3])
            elif find(l, 'Estimated Free Energy of Binding')>-1:
                d['binding_energy'] = float((ll)[7])
                if version>=4.0: d['energy'] = float((ll)[7])
            elif find(l, 'vdW + Hbond + desolv Energy')>-1:
                d['vdw_hb_desolv_energy'] = float((ll)[8])
            elif find(l, 'Electrostatic Energy')>-1:
                if (ll[1] == 'Electrostatic'):
                    d['electrostatic_energy'] = float((ll)[4])
                elif (ll[1] == 'Intermol.'):
                    d['electrostatic_energy'] = float((ll)[5])
            elif find(l, 'Moving Ligand-Fixed Receptor')>-1:
                d['moving_ligand_fixed_receptor'] = float((ll)[5])
            elif find(l, 'Moving Ligand-Moving Receptor')>-1:
                d['moving_ligand_moving_receptor'] = float((ll)[5])
            elif find(l, 'Total Internal Energy')>-1:
                d['total_internal'] = float((ll)[7])
            elif find(l, 'Internal Energy Ligand')>-1:
                d['ligand_internal'] = float((ll)[5])
            elif find(l, 'Internal Energy Receptor')>-1:
                d['receptor_internal'] = float((ll)[5])
            elif find(l, 'Torsional Free Energy')>-1:
                d['torsional_energy'] = float((ll)[6])
            elif find(l, 'Estimated Inhibition Constant')>-1:
                if find(l, 'N/A') < 0:
                    d['inhib_constant'] = float((ll)[6])
                    if version>=4.0:
                        d['inhib_constant_units'] = ll[7]
            elif find(l, 'Final Docked Energy')>-1:
                d['docking_energy'] = float((ll)[5])
            elif find(l, 'Final Intermolecular Energy')>-1:
                d['intermol_energy'] = float((ll)[6])
            elif find(l, 'Final Internal Energy of Ligand')>-1:
                d['internal_energy'] = float((ll)[8])
            elif find(l, 'Final Internal Energy')>-1:
                d['internal_energy'] = float((ll)[6])
            elif find(l, 'Torsional Free Energy')>-1:
                d['torsional_energy'] = float((ll)[6])
            elif find(l, 'NEWDPF tran0')>-1:
                d['trn_x'] = float(ll[3])
                d['trn_y'] = float(ll[4])
                d['trn_z'] = float(ll[5])
            elif find(l, 'NEWDPF quat0')>-1:
                d['qtn_nx'] = float(ll[3])
                d['qtn_ny'] = float(ll[4])
                d['qtn_nz'] = float(ll[5])
                d['qtn_ang_deg'] = float(ll[6])
            elif find(l, 'NEWDPF quaternion0')>-1:
                d['quaternion_nx'] = float(ll[3])
                d['quaternion_ny'] = float(ll[4])
                d['quaternion_ny'] = float(ll[5])
                d['quaternion_nw'] = float(ll[6])
            elif find(l, 'NEWDPF dihe0')>-1:
                angList = []
                for n in ll[3:]:
                    angList.append(float(n))
                d['torsion_values'] = angList
                d['num_torsions'] = len(angList)
            elif find(l, 'Intermol. vdW + Hbond Energy ')>-1:
                #AD4 specific model information:
                # USER         Intermol. vdW + Hbond Energy   =  -14.63 kcal/mol
                d['vdw_energy'] = float((ll)[7])
            elif find(l, 'Intermol. Electrostatic Energy')>-1:
                #USER         Intermol. Electrostatic Energy =   -0.62 kcal/mol
                d['estat_energy'] = float((ll)[5])
            elif find(l, '(3) Torsional Free Energy')>-1:
                #USER    (3) Torsional Free Energy           =   +3.84 kcal/mol
                d['torsional_energy'] = float((ll)[6])
            elif find(l, "(4) Unbound System's Energy")>-1:
                #USER    (4) Unbound System's Energy         =   -0.85 kcal/mol
                try:
                    d['unbound_energy'] = float((ll)[7])
                except:
                    d['unbound_energy'] = float((ll)[6])
            elif find(l, 'ATOM')>-1:
                coords.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                try:
                    vdW.append(float(l[54:60]))
                except:
                    #print 'vdw:', l[54:60], ' RAISED!'
                    vdW.append(0.0)
                try:
                    Elec.append(float(l[60:66]))
                except:
                    #print 'estat:', l[60:66], ' RAISED!'
                    Elec.append(0.0)
                if len(l)>77:
                    try:
                        binding_energy2 = float(l[70:76])
                    except:
                        pass
            elif find(l, 'HETA')>-1:
                coords.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                try:
                    vdW.append(float(l[54:60]))
                except:
                    vdW.append(0.0)
                try:
                    Elec.append(float(l[60:66]))
                except:
                    Elec.append(0.0)
            d['coords'] = coords
            d['vdw_energies'] = vdW
            d['estat_energies'] = Elec
            d['total_energies'] = Numeric.array(Numeric.array(vdW)+Numeric.array(Elec)).tolist()
            if binding_energy2 and not d.has_key('binding_energy'):
                d['binding_energy'] = binding_energy2
                d['docking_energy'] = binding_energy2
        return d



