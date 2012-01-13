#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2002
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/TrajParser.py,v 1.2 2003/08/29 17:55:18 sophiec Exp $
#
# $Id: TrajParser.py,v 1.2 2003/08/29 17:55:18 sophiec Exp $
#
#
#
#
#
#
#

"""
This Object parses the result of a trajectory AutoDock run and returns a dictionary. 

"""
import os
from string import find, join, replace, split, rfind
import re
from MolKit.pdbParser import PdbqParser

from AutoDockTools.ResultParser import ResultParser


class TrajParser(ResultParser):
    """ reads log from a AutoDock traj docking and return structured data"""

    keywords = ResultParser.keywords + [
        'run',
        'cycle',
        'temp',
        'acc_rej_code',
        #'intermol_energy',  #(1) 
        #'internal_energy',  #(2) NB: 1+2->final docked energy
        #'e_total',          #1+2 THIS IS THE SAME AS docking_energy????

        ]
        

    def __init__(self, dlgFile=None):
        """selected dlgFile,ok sets which docked conformations to show"""
        ResultParser.__init__(self)
        self.filename = dlgFile
        #set up dict here
        if dlgFile:
            self.filename = os.path.basename(dlgFile)
            self.parse(dlgFile)


    def parse(self, filename):
        """
        to parse uses keys:
            'ntorsions',
            'run',
            'cycle',
            'temp',
            'state'
        after parsing: 
            self.clist is list of dictionaries for states,
        """
        self.filename = filename
        #reset
        self.clist = []
        dlgptr = open(filename, 'r')
        allLines =  dlgptr.readlines()
        #first 4 lines are ntorsions, run, cycle and temp
        for i in range(len(allLines)):
            l = allLines[i]
            ll = split(l)
            if 'ntorsions'==ll[0]:
                self.num_torsions = int(ll[1])
            elif 'run' == ll[0]:
                self.run = int(ll[1])
            elif 'cycle' == ll[0]:
                self.cycle = int(ll[1])
            elif 'temp' == ll[0]:
                self.temp = float(ll[1])
            elif 'state' == ll[0]:
                self.getTrajState(allLines[i:])


    def getTrajState(self, lines):
        # BUILD A DICTIONARY and put it in clist
        d = {}
        d['num_torsions'] = self.num_torsions
        #lines[i] has format:
        #state 1 A e_total x y z qx qy qz qw
        #torsion 1
        # etc
        #WHEREAS in 'test-1.dlg' State= + 17 items: 3 trans, 4quat + ndihe(10) torsions
        xx = split(lines[0])

        # remove possible punctuation
        for ind in range(len(xx)):
            if xx[ind][-1]==',': xx[ind] = xx[ind][:-1]
            if xx[ind][-1]=='.': xx[ind] = xx[ind][:-1]

        d['cycle'] = int(xx[1])
        d['acc_rej_code'] = xx[2]
        d['docking_energy'] = float(xx[3])
        d['trn_x'] = float(xx[-7])
        d['trn_y'] = float(xx[-6])
        d['trn_z'] = float(xx[-5])
        d['qtn_nx'] = float(xx[-4])
        d['qtn_ny'] = float(xx[-3])
        d['qtn_nz'] = float(xx[-2])
        d['qtn_ang_deg'] = float(xx[-1])

        angList = []
        #NB: here torsions are NOT in the same line
        for i in range(1, self.num_torsions+1):
            angList.append(float(split(lines[i])[0]))
        d['torsion_values'] = angList
        self.clist.append(d)
        return 

