#########################################################################
#
# Date: Nov 2001 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################

from NetworkEditor.items import NetworkNode

from AutoDockTools.filterLigands import FilterLigands
from AutoDockTools.VisionInterface.Adt.LigandDB import LigandDB

class FilterLigandsNode(NetworkNode):
    """
    A node that takes a ligandDB, applies filtering and output a LigandList
    """

    attr = [
        'HbAMin', 'HbAMax', 'HbDMin', 'HbDMax', 'MWMin', 'MWMax',
        'NatMin', 'NatMax', 'TORSDOFMin', 'TORSDOFMax' ]

    def setParamPanelWidget(self, lfilter):
        """
        Set the widgets of the parameter panel to the filtering values without
        triggering network execution
        """
        for i in range(10):
            w = self.inputPorts[2+i].widget
            value = getattr(lfilter, self.attr[i])
            if value in [None, 'None']:
                w.configure(showLabel=0)
            else:
                w.configure(showLabel=1)
                w.set(value, run=0)


    def __init__(self, name='filterLigands', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        # create filter object
        lfilter = self.lfilter = FilterLigands()
        lfilter.setFilter('default')
        
        ip = self.inputPortsDescr
        ip.append(datatype='LigandDB', name='ligands')
        ip.append(datatype='string', name='filterMode')

        self.widgetDescr['filterMode'] = {
            'class':'NEComboBox', 'master':'node',
            'choices':['None', 'default', 'Lipinski-like', 'Drug-like', 'Drug-like frag', 'Custom'],
            'fixedChoices':True,
            'initialValue':'default',
            'entryfield_entry_width':8,
            'labelGridCfg':{'sticky':'w'},
            'widgetGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'filter type:'}}

        labels = [
            'Min num. of h-bond donors:',
            'Max num. of h-bond donors:',
            'Min num. of h-bond acceptors',
            'Max num. of h-bond acceptors',
            'Min molecular weight',
            'Max molecular weight',
            'Min num. of heavy atoms',
            'Max num. of heavy atoms',
            'Min num. of rotatable bonds',
            'Max num. of rotatable bonds'
            ]

        for i, param in enumerate( ['hbd_min', 'hbd_max', 'hba_min', 'hba_max',
                                    'mw_min', 'mw_max', 'nat_min', 'nat_max',
                                    'torsdof_min', 'torsdof_max'] ):
            ip.append(datatype='int', name=param)
            value = getattr(lfilter, self.attr[i])
            
            self.widgetDescr[param] = {
                'class':'NEThumbWheel','master':'ParamPanel', 
                'width':75, 'height':12, 'oneTurn':10, 'type':'int',
                'wheelPad':0, 'initialValue':value, 'min':0,
                'labelCfg':{'text':labels[i]} }

        op = self.outputPortsDescr
        op.append(datatype='LigandDB', name='ligands')
        op.append(datatype='int', name='num_accepted')
        op.append(datatype='int', name='num_rejected')

        code = """def doit(self, ligands, filterMode='default', hbd_min=None, hbd_max=None, hba_min=None, hba_max=None, mw_min=None, mw_max=None, nat_min=None, nat_max=None, torsdof_min=None, torsdof_max=None):

    lfilter = self.lfilter

    if filterMode != "Custom":
        lfilter.setFilter(filterMode)
    else:
        kwargs = {"HbDMin": hbd_min,
                  "HbDMax": hbd_max,
                  "HbAMin": hba_min,
                  "HbAMax": hba_max,
                  "MWMin": mw_min,
                  "MWMax": mw_max,
                  "NatMin": nat_min,
                  "NatMax": nat_max,
                  "TORSDOFMin": torsdof_min,
                  "TORSDOFMax": torsdof_max}
        lfilter.setFilterRanges('Custom', **kwargs)

    if self.inputPorts[1].hasNewValidData(): # new mode -> update widgets
        self.setParamPanelWidget(lfilter)

    print "Running Filter node with parameters:"
    print "min number of hydrogen bond donors: " + str(lfilter.HbDMin)
    print "max number of hydrogen bond donors: " + str(lfilter.HbDMax)
    print "min number of hydrogen bond acceptors: " + str(lfilter.HbAMin)
    print "max number of hydrogen bond acceptors: " + str(lfilter.HbAMax)
    print "min molecular weight: " + str(lfilter.MWMin)
    print "max molecular weight: " + str(lfilter.MWMax)
    print "min number of heavy atoms: " + str(lfilter.NatMin)
    print "max number of heavy atoms: " + str(lfilter.NatMax)
    print "min number of rotatable bonds: " + str(lfilter.TORSDOFMin)
    print "max number of rotatable bonds: " + str(lfilter.TORSDOFMax)

#    accepted, rejected = lfilter.filterTable(ligands.propertyTable)

    accepted, rejected = lfilter.filterTable(ligands.propertyTable, ligands.accepted)

    ligands.SetAcceptedLigands(accepted.filenames)

    if len(accepted.filenames) > 2500:
        print "ERROR: Sorry, we cannot send your virtual screening job to the server because the number of accepted ligands from the filter is greater than the maximum number of ligands (2500) the server can support"
        return 'stop'

    self.outputData(ligands=ligands, num_accepted = len(accepted.filenames), num_rejected=len(rejected))
"""
        self.setFunction(code)
