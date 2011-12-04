#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_clustering_histogram_postscript.py,v 1.3 2007/10/08 18:14:58 rhuey Exp $
#
# $Id: write_clustering_histogram_postscript.py,v 1.3 2007/10/08 18:14:58 rhuey Exp $
#
import os, glob, Tkinter, numpy.oldnumeric as Numeric
from MolKit import Read

from AutoDockTools.Docking import Docking
from AutoDockTools.interactiveHistogramGraph import InteractiveHistogramGraph
from AutoDockTools.histogram import HistogramRI
from mglutil.math.rmsd import RMSDCalculator




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: write_clustering_histogram_postscript.py -d directory"
        print
        print "    Description of command..."
        print "         -d     directory"
        print "    Optional parameters:"
        print "        [-t]    rmsd tolerance (default is 1.0)"
        print "        [-o]    output filename"
        print "                      (default is 'directory.ps')"
        print "        [-v]    verbose output"
        print "                      (default is leave active)"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:o:t:vh')
    except getopt.GetoptError, msg:
        print 'write_clustering_histogram_postscript.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None
    #-t: rms_tolerance
    rms_tolerance =  1.0
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = ""


    #'d:o:t:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            directory = a
            if verbose: print 'set directory to ', a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print 'set rms_tolerance to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    print 'directory=', directory
    if directory is None:
        print 'write_clustering_histogram_postscript: directory must be specified.'
        usage()
        sys.exit()
    if outputfilename=="":
        outputfilename = directory + '_' + str(rms_tolerance)+ '.ps'

    #read all the docking logs in as one Docking
    dlg_list = glob.glob(directory + '/*.dlg')
    d = Docking()
    for dlg in dlg_list:
        d.readDlg(dlg)
        if len(d.ch.conformations)>=50:
            print 'stop reading dlgs at', len(d.ch.conformations)
            break
    mol = d.ligMol
    crds = mol.allAtoms.coords[:]
    d.clusterer.rmsTool = RMSDCalculator(crds)
    d.clusterer.rmsToolRef = '0'   # for ADT only?
    #d.clusterer_dict['ats_binding'] = d.clusterer #for ADT only??
    d.clusterer.make_clustering(rms_tolerance) 
    s = d.ch.conformations
    elist = []
    for c in s: 
        if c.docking_energy:
            elist.append(c.docking_energy)
        elif c.binding_energy: 
            elist.append(c.binding_energy)
        elif c.energy: 
            elist.append(c.energy)
    if not len(elist): 
        print 'No energies available for docking.ch.conformations'
        exit()
    r = Tkinter.Tk()
    dataList = []
    reverseList = []
    rLctr = 0
    confL = d.ch.conformations
    e = d.clusterer.energy_used
    #for l in mol.cluSEQ:
    for l in d.clusterer.clustering_dict[rms_tolerance]:
        dataList.append([l[0].energy, len(l)])
        reverseList.append(range(rLctr, rLctr+len(l)))
    mol.elist = Numeric.array(elist)
    mol.r = [Numeric.minimum.reduce(mol.elist), 
                Numeric.maximum.reduce(mol.elist)]
    mol.nbins = Tkinter.IntVar()
    mol.nbins.set(10)
    mol.min = Tkinter.StringVar()
    mol.min.set(str(mol.r[0]))
    mol.max = Tkinter.StringVar()
    mol.max.set(str(mol.r[1]))

    r = (float(mol.min.get()), float(mol.max.get()))
    mol.ehist = HistogramRI(mol.elist,mol.nbins.get(),range=r)
    mol.ehist.createReverseIndex()
    nodeList = mol.ehist.array
    tstr = mol.name + ' histogram'
    top = Tkinter.Toplevel()
    top.title(tstr)
    mol.ehist
    #top = Tkinter.Toplevel()
    xlabel = 'ENERGY'
    mol.clustNB = InteractiveHistogramGraph(mol.name,
            master=top, nodeList = dataList, reverseIndex=reverseList,
            label_text=mol.name + ':' + str(rms_tolerance) + ' rms', xlabel_text=xlabel, 
            ylabel_text='#\nC\nO\nN\nF\nO\nR\nM\nA\nT\nI\nO\nN\nS')
    mol.clustNB.draw.update()
    mol.clustNB.draw.postscript({'file':outputfilename, 'colormode':'color'})
    top.update_idletasks()

# To execute this command type:
# write_clustering_histogram_postscript.py -d directory -t rmsd tolerance -o outputfilename

