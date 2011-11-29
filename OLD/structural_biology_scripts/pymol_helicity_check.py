# pymol_helicity_check.py
# Copyright (c) 2006-2007 Julien Lefeuvre <lefeuvrejulien@yahoo.fr>
#
 
"""
Pymol plugin for checking helicity type
 
helicity_check() takes as input a selection ('sele' by default)
of at least 5 amino acids and computes the distances between
O(i) - N(i+3)
O(i) - N(i+4)
O(i) - N(i+5)
See for further info:
Protein Sci ROHL and DOIG 5 (8) 1687
'Models for the 3(10)-helix/coil, pi-helix/coil,
and alpha-helix/3(10)-helix/coil transitions in isolated peptides.'
 
uses:
*in the pymol console:
  >run pymol_helicity_check.py
    ----> select some consecutive amino acids
           - this is nicely done with the Display->Sequence tool
  >helicity_check()
*installing helicity_check
  copy pymol_helicity_check.py in $PYMOL_INSTALL_DIR/modules/pmg_tk/startup
  launch Pymol: you now have a new option in the Plugin menu
 
helicity_check uses gnuplot (http://www.gnuplot.info) to display its results
As a consequence gnuplot needs to be installed.
 
This plugin was tested on linux only, it my need some modifications to run on
other OSes (hints: launching gnuplot and path to dumpfile)
"""
 
__author__ =    "Julien Lefeuvre <lefeuvrejulien@yahoo.fr>"
__version__ =   "1.0"
__date__ =      "2007-04-02"
__copyright__ = "Copyright (c) 2007 %s. All rights reserved." % __author__
__licence__ =   "BSD"
 
from pymol import cmd
from math import sqrt
import sys
import os
import subprocess
import time
 
def __init__(self):
    """init function in order to have a nice menu option in Pymol"""
    self.menuBar.addmenuitem('Plugin', 'command', 'Helicity Check',
             label='Helicity Check', command = lambda: helicity_check())
 
 
class Residue(object):
 
    def __init__(self):
        self.name=None
        self.index=None
        self.Ocoord=None
        self.Ncoord=None
 
 
def calc_distON(Ocoord,Ncoord):
    """return the distance between 2 atoms given their coordinates"""
    sum = 0
    for o, n in zip(Ocoord, Ncoord):
        sum += (o - n)**2
    return sqrt(sum)
 
 
def helicity_check(selection='sele'):
    """calcultate distance O[res i]-N[res i+3]
                           O[res i]-N[res i+4]
                           O[res i]-N[res i+5]
    """
    seq_model = cmd.get_model(selection) #get info from selection
    res_lim = seq_model.get_residues()
 
    if len(res_lim)<5:
        sys.stderr.write("\nPlease select at least 5 residues\n")
        return
 
    atom_list = seq_model.atom
    res_data=[]
 
    for start,end in res_lim:   #extract the data we are interested in
        res=Residue()
        for atom in atom_list[start:end]:
            if atom.name == 'N':
                res.name = atom.resn
                res.index = atom.resi
                res.Ncoord = atom.coord
            elif atom.name == 'O':
                res.Ocoord = atom.coord
        if res.Ocoord and res.Ncoord and res.name and res.index:
            res_data.append(res)
        else:
            sys.stderr.write("\nPlease select complete protein residues\n")
            return
 
    res_list = [int(res.index) for res in res_data]
 
    if res_list != range(res_list[0], res_list[-1]+1):
        sys.stderr.write("\nPlease select a unbrocken residue sequence\n")
        return
 
    distON3 = []
    distON4 = []
    distON5 = []
    distONs = [distON3, distON4, distON5]
 
    for i,res in enumerate(res_data[:-5]): #distances calculations
        resis = res_data[i+3:i+6]
        for resi, distONi in zip(resis, distONs):
            distONi.append(calc_distON(res.Ocoord, resi.Ncoord))
 
    dump = os.tmpnam()+'.dat'
    dumpfile = file(dump, 'w')
 
    sys.stdout.write('\n#Distances O(i)---N(i+n)\n'
           '#ResNum , d(O(i)-N(i+3)) , d(O(i)-N(i+4)) , d(O(i)-N(i+4))\n')
    for i, d3, d4, d5 in zip(res_list, distON3, distON4, distON5):
        #writing console output
        sys.stdout.write(
              '  %i ,      %f ,       %f ,       %f \n'%(i, d3, d4, d5))
        #writing data to a dump file for use by gnuplot
        dumpfile.write(
              '  %i       %f        %f        %f \n'%(i, d3, d4, d5))
    dumpfile.flush()
 
    #launch a gnuplot window to show the distances
    gnuplotcmd = subprocess.Popen(['/usr/bin/gnuplot'], shell=True,
                               stdin=subprocess.PIPE)
    gnuplotcmd.stdin.write('set autoscale\n')
    gnuplotcmd.stdin.write("plot "
         "'%s' using 1:2 title 'd(O(i)-N(i+3))' with lines, "
         "'%s' using 1:3 title 'd(O(i)-N(i+4))' with lines, "
         "'%s' using 1:4 title 'd(O(i)-N(i+5))' with lines\n'"
                          % (dump, dump, dump))
    time.sleep(3)
    dumpfile.close()
    os.remove(dump)
