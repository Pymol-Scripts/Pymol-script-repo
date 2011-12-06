#
# Last modified on Mon Mar  4 14:35:36 PST 2002 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/ResultParser.py,v 1.3 2009/01/09 19:53:49 rhuey Exp $
#




class ResultParser:
    """Base class of result parsers. Supply your own parseline
    function which appends a conformation dictionary to self.clist
    """
    def __init__(self):
        self.clist = []

    # add more keywords in subclass if you like 
    keywords = [
        'cluster',            # number of cluster
        'cluster_rank',       # rank within cluster
        'rmsd_ref',           # distance to reference structure
        'rmsd_seed',          # distance to lowest energy conf. in cluster
        'binding_energy',     # estimated free energy of binding
        'docking_energy',     # final docked energy
        'internal_energy',
        'trn_x',              # translation x, y, z
        'trn_y',
        'trn_z',
        'qtn_nx',             # quaternion unit vector x, y, z
        'qtn_ny',
        'qtn_nz',
        'qtn_ang_deg',        # quaternion rotation angle
        'num_torsions',
        'torsion_values',
        'rseed1',             # the random number seeds for this conformation
        'rseed2',
        ]

    def parse(self, filename):
        """
        """
        file_ptr = open(filename)
        
        self.clist = []
        for line in file_ptr.readlines():
            self.parseline(line)
        file_ptr.close()

        self.filename = filename
        
        return self.clist


    def parseline(self, line):
        """over ride me"""
        pass
    
