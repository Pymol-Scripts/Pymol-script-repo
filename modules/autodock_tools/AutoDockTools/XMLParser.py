#############################################################################
#
# Author: Ruth HUEY, William Lindstrom
#
# Copyright: M. Sanner TSRI 2005
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/XMLParser.py,v 1.14 2008/09/02 22:30:59 gillet Exp $
#
# $Id: XMLParser.py,v 1.14 2008/09/02 22:30:59 gillet Exp $
#
#
#
#
#
#
#

"""
This Object parses the xml result of an AutoDock operation. It builds a dictionary. 

"""
import os
from string import find, join, replace, split, rfind
import re

from AutoDockTools.ResultParser import ResultParser


class XMLParser(ResultParser):
    """ reads log from a AutoDock docking and return structured data"""
    
    keywords = ResultParser.keywords + [
        #'seed',   #rseed1, rseed2
        'dpf',
        #'free_NRG_binding',  #binding_energy
        'Ki',
        'Temp',
        'final_intermol_NRG', 
        #'internal_ligand_NRG', #internal_enrgy
        'torsional_free_NRG',
        'move',
        'about',
        #'tran0',  #trn_x, trn_y, trn_z
        #'quat0',  #qtn_nx, qtn_ny, qtn_nz, qtn_ang_deg
        #'ndihe',  #num_torsions
        #'dihe0',  # torsion_values
    ]


    def __init__(self, dlgFile=None, dpfFile=None):
        """selected dlgFile,ok sets which docked conformations to show"""
        ResultParser.__init__(self)
        self.filename = dlgFile
        self.version = 1.0
        if dlgFile:
            self.filename = os.path.basename(dlgFile)
            self.parse(dlgFile)
        if dpfFile:
            self.dpf = dpfFile
            

    def parse(self, filename):
        """
        uses key '<autodock>' to start matching:
        next uses '<runs>' to start capturing individual docked results 
        finally captures '</autodock>' to end
        after parsing: 
        """
        self.filename = filename
        #reset
        dlgptr = open(filename, 'r')
        allLines = self.allLines = dlgptr.readlines()
        self.clusterRecord = None
        #print "calling match with ", len(allLines)
        self.match(allLines)


    def getReDict(self):
        if hasattr(self, 'reDict'):
            for k, d in self.reDict.items():
                d['lines'] = []
            return
        self.reDict = {}
        self.reKeys = [
        '\t<version>',
        '\t<autogrid_version>',
        '\t<output_xml_version>',
        '\t<run_requested>',
        '\t<runs>',
            ]
        self.reFuncs = [ 
            self.set_AD_version,
            self.set_AG_version,
            self.set_XML_version,
            self.set_runs_requested,
            self.get_runs,
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
        #for i in range(5):
        for i in range(len(allLines)):
            item = allLines[i]
            #print "item=", item
            #if find item, mark it found + don't test
            for k in self.reKeys:
                d = self.reDict[k]
                m = d['re'].match(item)
                if m:
                    #print "matched ", k
                    d['lines'].append(item) 
                    break
        for k in self.reKeys:
            d = self.reDict[k]
            lines = d['lines']
            apply(d['func'], (lines,), {})


    def set_AD_version(self,lines):
        if len(lines):
            for l in lines:
                if find(l, '<version>')>-1:
                    ll = l.split('>')
                    if len(ll)>0:
                        lll = ll[1].split('<')
                        self.version = float(lll[0])
                    else:
                        print "problem autodock version found!"
                        self.version = 4.03
                        break
                    #print "ad version=", self.version
                    break
        else:
            print "no autodock version found!"
            self.version = 4.03


    def set_AG_version(self,lines):
        if len(lines):
            for l in lines:
                if find(l, '<autogrid_version>')>-1:
                    ll = l.split('>')
                    if len(ll)>0:
                        lll = ll[1].split('<')
                        self.autogrid_version = float(lll[0])
                    else:
                        print "problem autogrid version found!"
                        self.autogrid_version = 4.03
                        break
                    #print "ag version=", self.autogrid_version
                    break
        else:
            print "no autogrid version found!"
            self.autogrid_version = 4.03


    def set_XML_version(self,lines):
        if len(lines):
            for l in lines:
                if find(l, '<output_xml_version>')>-1:
                    ll = l.split('>')
                    if len(ll)>0:
                        lll = ll[1].split('<')
                        self.xml_version = float(lll[0])
                    else:
                        print "problem xml version version found!"
                        self.xml_version = 0.10
                        break
                    #print "xml version=", self.xml_version
                    break
        else:
            print "no xml version found!"
            self.xml_version = 0.10


    def set_runs_requested(self, lines):
        if len(lines):
            for l in lines:
                if find(l, '<run_requested>')>-1:
                    ll = l.split('>')
                    if len(ll)>0:
                        lll = ll[1].split('<')
                        self.run_requested = int(lll[0])
                    else:
                        print "problem with run requested found!"
                        self.run_requested = 1
                        break
                    #print "run requested =", self.run_requested
                    break


    def get_runs(self, lines):
        #print "in get runs with lines=", lines
        if len(lines):
            if lines[0].find('<runs>')>-1:
                ind = self.allLines.index(lines[0])
                run_lines = []
                for l in self.allLines[ind:]:
                    if l.find('</runs>')>-1:
                        return
                    elif l.find('</run>')>-1:
                        #ends with </run>
                        #print "end of a run!"
                        self.process_run(run_lines)
                        run_lines = []
                    else:
                        #starts with <run id="   1">
                        #accummulate lines for each run
                        #process these lines and 
                        run_lines.append(l)


    def get_floats(self, line):
        #print "in get_floats with ", line
        ll = line.split('>')[1]
        lll = ll.split('<')
        return map(float, lll[0].split())


    def get_ints(self, line):
        #print "in get_ints with ", line
        ll = line.split('>')[1]
        lll = ll.split('<')
        return map(int, lll[0].split())
        

    def process_run(self, run_lines):
        #seed, dpf, free_NRG_binding, Ki, Temp, final_intermol_NRG
        #internal_ligand_NRG, torsional_free_NRG, move, about, tran0, quat0,
        #ndihe, dihe0
        #print "in process_run with ", run_lines
        dihe_list = [] # initial for rigid ligands (no ndihe tag)
        for l in run_lines:
            if l.find('<run id="')>-1:
                id = int(l.split('"')[1])
            elif l.find('<seed>')>-1:
                seed1, seed2 = self.get_ints(l)
            elif l.find('<dpf>')>-1:
                #print l
                #print l[7:-7]
                dpf = l[7:-7]
                #IGNORE dpf in file
                if hasattr(self, 'dpf'):
                    #print "ignoring dpf in file"
                    dpf = self.dpf
            elif l.find('<free_NRG_binding>')>-1:
                free_NRG_binding = self.get_floats(l)[0]
            elif l.find('<Ki>')>-1:
                Ki = self.get_floats(l)[0]
            elif l.find('Temp')>-1:
                Temp = self.get_floats(l)[0]
            elif l.find('final_intermol_NRG')>-1:
                final_intermol_NRG = self.get_floats(l)[0]
            elif l.find('internal_ligand_NRG')>-1:
                internal_ligand_NRG = self.get_floats(l)[0]
            elif l.find('torsonial_free_NRG')>-1:
                torsional_free_NRG = self.get_floats(l)[0]
            elif l.find('move')>-1:
                self.ligand = l.split('>')[1].split('<')[0]
            elif l.find('about')>-1:
                about = self.get_floats(l)
            elif l.find('tran0')>-1:
                trans = self.get_floats(l)
            elif l.find('quat0')>-1:
                axisangle = self.get_floats(l)
            elif l.find('quaternion0')> -1:
                quaternion0 = self.get_floats(l)
            elif l.find('ndihe')>-1:
                ndihe = self.get_ints(l)[0]
                ind = run_lines.index(l)
                l = run_lines[ind+1]
                dihe_list = self.get_floats(l)
                rest_of_run_lines = run_lines[ind+2:]
                for l in rest_of_run_lines:
                    if l.find('</dihe0>')==-1:
                        more_dihe = self.int_floats(l)
                        dihe_list.extend(more_dihe)
                    else:
                        break
        d = {}
        d['id'] = id
        d['rseed1'] = seed1
        d['rseed2'] = seed2
        d['org_x'] = about[0]
        d['org_y'] = about[1]
        d['org_z'] = about[2]
        d['trn_x'] = trans[0]
        d['trn_y'] = trans[1]
        d['trn_z'] = trans[2]
        d['qtn_nx'] = axisangle[0]
        d['qtn_ny'] = axisangle[1]
        d['qtn_nz'] = axisangle[2]
        d['qtn_ang_deg'] = axisangle[3]
        d['num_torsions'] = d['ndihe'] = len(dihe_list)
        d['torsion_values'] = d['dihe0'] = dihe_list
        #XML specific values
        try:
            d['Ki'] = Ki
        except:
            print 'except on Ki, for run id=', id
            d['Ki'] = 0.0
        try:
            d['Temp'] = Temp
        except:
            print 'except on Temp, for run id=', id
            d['Temp'] = 0.0
        #intermol+internal+torsional
        d['binding_energy'] = d['free_NRG_binding'] = free_NRG_binding 
        d['intermol_energy'] = d['final_intermol_NRG'] = final_intermol_NRG
        d['internal_energy'] = d['internal_ligand_NRG'] = internal_ligand_NRG
        d['torsional_energy'] = d['torsional_free_NRG'] = torsional_free_NRG
        #print "d['torsional_energy']=", torsional_free_NRG
        d['quat0'] = axisangle
        d['quaternion0']= quaternion0
        d['tran0'] = trans
        d['about'] = about
        d['dpf'] = dpf
        d['xml'] = self.filename
        d['dlg'] = self.filename.replace('.xml','.dlg')
        
        if hasattr(self, 'dpf'):
            if dpf!=self.dpf:
                print "dpf mismatch"
                print dpf, ' vs ', self.dpf
            ##assert dpf == self.dpf
        else:
            self.dpf = dpf
        #print "appending new d"
        self.clist.append(d)
        #print "len(clist)=", len(self.clist)

