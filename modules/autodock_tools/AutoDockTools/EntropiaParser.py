#
# Last modified on Thu Mar 11 10:50:58 PST 2004 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/EntropiaParser.py,v 1.7 2004/03/11 19:17:04 lindy Exp $
#
#


import string
import types
from AutoDockTools.ResultParser import ResultParser


class EntropiaParser(ResultParser):
    """<class docstring>
    """
    # keywords are the keys of the dictionary used by the
    # ConformationHandler to instantiate a Conformation.
    keywords = ResultParser.keywords + [
        # add keywords here
        ]


    def __init__(self):
        ResultParser.__init__(self)


    def parseline(self, line):
        """Parse the given line and return append a dictionary
        onto self.clist, the superclass-declared list of conformations.
        """
        line_list = string.split(line)
        # make a heuristic decision about which keys to use!!
        if line_list[2] == repr(1.01) or line_list[2][:3] == repr(1.0):
            # version 1.01 Entropia result file
            dict = self._parseStateLineList(line_list, self.result_keys_v101)
        elif float(line_list[9]) == 1.0:
            # version 1.0 Entropia result file
            dict = self._parseStateLineList(line_list)
        else:
            print "Unparsable result file line: %s" % (line)
            return # without doing anything
        self.clist.append(dict)


    result_keys = [        # original version 1.00
        'output_id',
        'data_run_id',
        'dpf_id',
        'creation_dtime',
        'last_update_dtime',
        'ei_version',         # Entropia Interface version
        'ag_version',         # AutoGrid version
        'ad_version',         # AutoDock version
        'run_rank',           # rank of this run in cluster
        'run_number',         # number of this run in dpf
        'cluster_rank',       # rank of this cluster
        'cluster_size',       # number of conformations in this cluster
        'run_size',           # number of runs specified in dpf
        'rseed1',
        'rseed2',
        'rmsd',               # RMSD from reference structure
        'binding_energy',     # estimated free energy of binding
        'docking_energy',     # final docked energy
        'trn_x',              # translation x, y, z
        'trn_y',
        'trn_z',
        'qtn_nx',             # quaternion unit vector x, y, z
        'qtn_ny',
        'qtn_nz',
        'qtn_ang_deg',        # quaternion rotation angle
        'num_torsions',
        'torsion_values']


    result_keys_v101 = [   # for the 1.01 version of the Entropia Interface
        'data_run_id',
        'dpf_id',
        'ei_version',         # Entropia Interface version
        'ag_version',         # AutoGrid version
        'ad_version',         # AutoDock version
        'run_rank',           # rank of this run in cluster
        'run_number',         # number of this run in dpf
        'cluster_rank',       # rank of this cluster
        'cluster_size',       # number of conformations in this cluster
        'run_size',           # number of runs specified in dpf
        'rseed1',
        'rseed2',
        'rmsd',               # RMSD from reference structure
        'binding_energy',     # estimated free energy of binding
        'docking_energy',     # final docked energy
        'trn_x',              # translation x, y, z
        'trn_y',
        'trn_z',
        'qtn_nx',             # quaternion unit vector x, y, z
        'qtn_ny',
        'qtn_nz',
        'qtn_ang_deg',        # quaternion rotation angle
        'num_torsions',
        'torsion_values']


    type_conv = {
        types.IntType    : lambda x: int(x),
        types.FloatType  : lambda x: float(x),
        types.StringType : lambda x: x,
        'PercentType'    : lambda x: float(x[:-1]),
        'TimeType'       : lambda x: x
        }


    field_defn = {
    #    key                : (type, whitespace delimited sub-fields)
        'output_id'         : (types.IntType, 1),
        'data_run_id'       : (types.IntType, 1),
        'dpf_id'            : (types.IntType, 1),
        'creation_dtime'    : ('TimeType', 3),
        'last_update_dtime' : ('TimeType', 3),
        'ei_version'        : (types.FloatType, 1),
        'ag_version'        : (types.FloatType, 1),
        'ad_version'        : (types.FloatType, 1),
        'run_rank'          : (types.IntType, 1),
        'run_number'        : (types.IntType, 1),
        'cluster_rank'      : (types.IntType, 1),
        'cluster_size'      : (types.IntType, 1),
        'run_size'          : (types.IntType, 1),
        'rseed1'            : (types.IntType, 1),
        'rseed2'            : (types.IntType, 1),
        'rmsd'              : (types.FloatType, 1),
        'binding_energy'    : (types.FloatType, 1),
        'docking_energy'    : (types.FloatType, 1),
        'trn_x'             : (types.FloatType, 1),
        'trn_y'             : (types.FloatType, 1),
        'trn_z'             : (types.FloatType, 1),
        'qtn_nx'            : (types.FloatType, 1),
        'qtn_ny'            : (types.FloatType, 1),
        'qtn_nz'            : (types.FloatType, 1),
        'qtn_ang_deg'       : (types.FloatType, 1),
        'num_torsions'      : (types.IntType, 1),
        'torsion_values'    : (types.FloatType, None)
        }



    def _parseStateLineList(self, lineList = None, keys=None):
        """Return a dictionary of values.
        """
        if not keys:
            keys = self.result_keys
        
        # initalize the stateDict with all fields
        stateDict = {}
        for field in self.result_keys:
            stateDict[field] = None

        # now run through the specified fields
        subFldIx = 0
        for field in keys:

            # most fields just get converted in the first case
            if self.field_defn[field][1] == 1:
                stateDict[field] = self.type_conv[
                    self.field_defn[field][0]](lineList[subFldIx])
                subFldIx = subFldIx + 1

            # the time fields have spaces in them; handled here
            elif self.field_defn[field][0] == 'TimeType':
                timeList = []
                for subIx in range( self.field_defn[field][1]):
                    # just append the strings to the timeList
                    timeList.append(lineList[subFldIx])
                    subFldIx = subFldIx + 1
                stateDict[field] = timeList

            # this case handles the variable number of torsions
            elif field == 'torsion_values':
                torsionList = []
                for subIx in range( stateDict['num_torsions']):
                    # loop through the torsion_values
                    torsionList.append( self.type_conv[
                        self.field_defn[field][0]](lineList[subFldIx]))
                    subFldIx = subFldIx + 1
                stateDict[field] = torsionList

        # make sure there are no left over fields
        assert len(lineList) == subFldIx
        return stateDict









