"""
Authors: Stefano Forli, M Sanner

April 2010

copyright TSRI

CalculateProperties can be used to process a directory and create a property
table for the ligands in that directory.

example:
    calcProp = CalculateProperties(verbose=2)
    propT = calcProp.processDirectory(mydir, '*.pdbqt')

    or for a folder hirearchy
    propT = calcProp.processDirectoryTree(root, extension='.pdbqt')

propT is a PropertyTable object that stores the directory where the files are
and has a table attribute containing the actual properties.

propT.saveAsTextFile(filename) saves the property table

A PropertyTable object can be recreated from a file:
   pt = PropertyTable(filename='lib1.prop')

or from a URL
   pt = PropertyTable(url='http://www.scripps.edu/~sanner/collab/lib1.prop')


The FilterLigands can be used to filter ligands based on properties.
It creates a LigandList object which contains a list of full names (i.e.
path and file name) for each seelcted ligand. using full names enables
combining molecules from multiple databases, This object provides several
pre-set ranges of values for the properties corresponding to 4 commonly used
filters.

example:
    lfilter = FilterLigands()
    lfilter.setFilter('default')
    accepted, rejected = lfilter.filterTable(pt)
    print 'default filter in %f second accepted %d rejected %d'%(
        lfilter.filteringTime, len(accepted), len(rejected))
    print lfilter.getAtomTypes(calcProp.propertyTable, accepted)

"""

class PropertyTable:
    """
    """

    def __init__(self, folder=None, table=None, url=None, filename=None):
        self.folder = folder     # folder containing the ligand files
        self.table = table       #
        self.url = url        # folder containing the ligand files
        self.filename = filename   # folder containing the ligand files
        if url is not None:
            self.fillFromURL(url)
        elif filename is not None:
            self.fillFromFile(filename)
        else:
            assert folder is not None and isinstance(table, dict)


    def saveAsTextFile(self, filename):
        """
        None <- CalculateProperties.saveAsTextFile(filename)

        saves the propertyTable in a file as a dictionary called propertyTable

        """ 
        TORSDOF = 0
        HbD = 1
        HbA = 2
        MW  = 3
        Nat = 4
        NotStdAT = 5
        Atypes = 6
        nameLength = max( map(len, self.table.keys()))
        f = open(filename, 'w')
        pad = max(0, nameLength-7)
        f.write('folder %s\n'%self.folder)
        f.write('# file '+' '*pad+' #TORS HbD HbA  MW   Nat nstd Types\n')
        fmt = "%"+str(nameLength)+"s %3d  %3d %3d %5.2f  %3d  %d"
        for key, v in self.table.items():
            f.write(fmt%(key, v[TORSDOF], v[HbD], v[HbA], v[MW],
                    v[Nat], v[NotStdAT]))
            for t in v[Atypes]:
                f.write(' %2s'%t)
            f.write("\n")
        f.close()


    def fillFromFile(self, filename):
        f = open(filename)
        lines = f.read().split('\n')
        f.close()
        self._parseLines(lines)


    def fillFromURL(self, url):
        import urllib
        
        fp = urllib.urlopen(url)
        lines = fp.read().split('\n')
        fp.close()
        self._parseLines(lines)


    def _parseLines(self, lines):
        propTable = {}
        folder = None
        for line in lines:
            if len(line)==0 or line[0]=='#':
                continue
            elif line[:6]=='folder':
                self.folder = line[6:].strip()
            else:
                w = line.split()
                propTable[w[0]] = [ int(w[1]), int(w[2]), int(w[3]),
                                    float(w[4]), int(w[5]), int(w[6]), w[7:] ]
        self.table = propTable

##     def saveAsPythonCode(self, filename):
##         """
##         None <- CalculateProperties.saveAsPythonCode(filename)

##         filename isa to have a .py extension
##         saves the propertyTable in a file as a dictionary called propertyTable
##         """
##         from os import path
##         assert path.splitext(filename)[1]=='.py'
##         f = open(filename, 'w')
##         f.write('propertyTable = ')
##         f.write(str(self.table))
##         f.close()


class CalculateProperties:
    """
    This function parses the lines from a pdbqt file provided in records
    and creates a property record for the this molecule.
    the property record contains the following information:
        TORSDOF: int - number of rotatable bonds
        atype: list of strings - list of tom types in the molecule
        HbD: int - number of hydrogen bond donors
        HbA: int - number of hydrogen bond acceptors
        MW: float - molecular weight
        Nat: int - number of atoms
        NotStdAT: bool - contains non standard atoms

    the property records are stored in the propertyTable using the filename as
    as key
    
    Usage:
        cp = CalculateProperties(verbose=0)
        cp.processDirectory(dirname, '*.pdbqt')

    """
    
    molecularWeight = {   # MW
        'H'      :   1   ,
        'HD'     :   1   ,
        'HS'     :   1   ,
        'C'      :   12  ,
        'A'      :   12  ,
        'N'      :   14  ,
        'NA'     :   14  ,
        'NS'     :   14  ,
        'OA'     :   16  ,
        'OS'     :   16  ,
        'F'      :   19  ,
        'Mg'     :   24  ,
        'MG'     :   24  ,
        'P'      :   31  ,
        'SA'     :   32  ,
        'S'      :   32  ,
        'Cl'     :   35.4,
        'CL'     :   35.4,
        'Ca'     :   40  ,
        'CA'     :   40  ,
        'Mn'     :   55  ,
        'MN'     :   55  ,
        'Fe'     :   56  ,
        'FE'     :   56  ,
        'Zn'     :   65.4,
        'ZN'     :   65.4,
        'Br'     :   80  ,
        'BR'     :   80  ,
        'I'      :  126  ,
        'e'      :    0  ,
        'd'      :    0   
        }


    def __init__(self, verbose=0):
        """
        Constructor

        if verbose > 1: processing time will be printed
        if verbose > 10: properties for each molecule willbe printed
        """
        self.verbose = verbose
        self.pocessingTime = 0.0
        self.directory = None
        
        
    def getRecords(self, filename):
        """
        records <- CalculateProperties.getRecords(filename)

        read the file and returns a list of lines in the file
        """
        file = open(filename, 'r')
        lines = file.readlines()
        file.close()
        return lines


    def processDirectoryTree(self, root, extension='.pdbqt'):
        """
        PropertyTable <- CalculateProperties.processDirectoryTree(root, extension='.pdbqt')
        find all files matching extension in the folder tree starting at root
        """
        import os, time
        
        splitext = os.path.splitext
        join = os.path.join
        t1 = time.time()
        propT = {}
        off = len(root)+1
        for rootl, dirs, files in os.walk(root):
            for f in files:
                if splitext(f)[1]==extension:
                    filename = join(rootl, f)
                    lines = self.getRecords(filename)
                    key = filename[off:] # exclude root
                    self.calcProp(lines, key, propT)
                    
        propTab = PropertyTable(folder=root, table=propT)

        self.processingTime = time.time()-t1
        if self.verbose>1:
            print 'Processed %d molecules in %.2f seconds'%(
                len(propT), self.processingTime)
        return propTab


    def processDirectory(self, dirname, selector='*'):
        """
        PropertyTable <- CalculateProperties.processDirectory(dirname, selector='*')
        """
        from glob import glob
        import os, time
        
        t1 = time.time()
        propT = {}
        filenames = glob(os.path.join(dirname, selector))
        for filename in filenames:
            lines = self.getRecords(filename)
            # the the filename without extension as the key
            key = os.path.splitext(os.path.basename(filename))[0]
            self.calcProp(lines, key, propT)

        propTab = PropertyTable(folder=dirname, table=propT)

        self.processingTime = time.time()-t1
        if self.verbose>1:
            print 'Processed %d molecules in %.2f seconds'%(
                len(propT), self.processingTime)
        return propTab
    

    def calcProp(self, records, key, propT):
        """
        """
        current_atypes = {}
        BAD_ATOM_TYPE = False
        MW  = 0
        HbD = 0
        HbA = 0
        Nat = 0
        status = True
        TORSDOF = 0

        hbd_h = []
        hbd_candidate = []

        # Calculate all the properties
        for line in records:
            if 'TORSDOF' in line: # not a defined location ??
                TORSDOF = int(line[8:])

            if line[0:6] == 'HETATM' or line[0:4] == 'ATOM':
                # Nat += 1 # Consider to remove Hydrogens? (search on PubMed)
                atype = line[76:-1].strip()
                if atype not in current_atypes:
                    current_atypes[atype] = True

                # Hb acceptor
                if atype == "OA" or atype == "NA" or atype == "SA":
                    HbA += 1

                # Hb donor preparation
                if atype == "HD":
                    #capture the hydrogens that could be bond to the Hb donor...
                    hbd_h.append( (float(line[30:38]), float(line[38:46]),
                                   float(line[46:54]) ) )
                else:
                    # count heavy atoms
                    Nat += 1

                if atype == "N" or atype == "O" or atype == "OA" or atype == "NA":
                    hbd_candidate.append( (float(line[30:38]),
                                   float(line[38:46]), float(line[46:54])) )
                try:
                    # add the atomic weight to the total MW
                    MW += self.molecularWeight[atype]
                except:
                    #MW += 10000 # check this if it's reasonable
                    MW += 0 # check this if it's reasonable
                    BAD_ATOM_TYPE = True
                    status = False # non-standard atom types are rejected by default

            # identify HBD by checking if N/O's are bound to H's
            dcut2 = 1.1 *1.1
            for atom in hbd_candidate:
                x1,y1,z1 = atom
                for hydrogen in hbd_h:
                    x2,y2,z2 = hydrogen
                    d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) 
                    if d2 <= dcut2:
                        HbD += 1
                        break

        propT[key] = [ TORSDOF, HbD, HbA, MW, Nat, BAD_ATOM_TYPE, 
                       current_atypes.keys() ]
#            {
#            "Atypes": current_atypes.keys(),
#            "TORSDOF": TORSDOF,
#            "HbD": HbD,
#            "HbA": HbA,
#            "MW": MW,
#            "Nat": Nat,
#            "NotStdAT":BAD_ATOM_TYPE
#            }

        if self.verbose > 10:
            print "tors %d hbd:%d hba: %d mw: %.2f nat:%d nstd: %d, types: %s\n"%(
                TORSDOF, HbD, HbA, MW, Nat, BAD_ATOM_TYPE,
                current_atypes.keys()) 
           

       


class FilterLigands:
    """
    Class for filterting a propertyTable

    usage:
        ligFilter = FilterLigands()
        ligFilter.setFilter('default')
        accepted, rejected = ligFilter.filterTable(propertyTable)

    setFilter() accepts the following filters:
        'default', 'Lipinski-like', 'Drug-like', 'Drug-like frag'
        It also can be passed 'None' or None to prevent filtering

    setFilterRanges(filterType, **kw) can be called to set filtering values
    to special ranges.The .filter attribute is set to 'Custom' in this case.
    """
    filterNames = ['default', 'Lipinski-like', 'Drug-like', 'Drug-like frag']

    def __init__(self):
        self.filter = None # stores the filter type
        self.setFilter('default')
        self.rejectNStdAt = True
        self.filteringTime = 0


    def setFilterRanges(self, filterType, **kw):

        if filterType in self.filterNames+['None', None]:
            self.filter=filterType
        else:
            self.filter = 'Custom'

        if len(kw)==0: # no options means reset values to default
            self.HbDMin = 0# HydBond DONORS
            self.HbDMax = 99
            self.HbAMin = 0 # HydBond ACCEPTOR
            self.HbAMax = 99
            self.MWMin = 0 # Molecular weight
            self.MWMax = 99999
            self.NatMin = 0 # Number of heavy atoms
            self.NatMax = 999
            self.TORSDOFMin = 0
            self.TORSDOFMax = 32
        else:
            for key, value in kw.items():
                if key=='HbDMin':
                    self.HbDMin = value
                elif key=='HbDMax':
                    self.HbDMax = value
                elif key=='HbAMin':
                    self.HbAMin = value
                elif key=='HbAMax':
                    self.HbAMax = value
                elif key=='MWMin':
                    self.MWMin = value
                elif key=='MWMax':
                    self.MWMax = value
                elif key=='NatMin':
                    self.NatMin = value
                elif key=='NatMax':
                    self.NatMax = value
                elif key=='TORSDOFMin':
                    self.TORSDOFMin = value
                elif key=='TORSDOFMax':
                    self.TORSDOFMax = value
                else:
                    raise ValueError("bad key %s for setting value range, expected: HHbDMin, HbDMax, HbAMin, HbAMax, MWMin, MWMax, NatMin, NatMax, TORSDOFMin, TORSDOFMax")


    def setFilter(self, mode='default'):
        """
        Configure the filtering ranges for various types of filters
        mode can be 'default', 'Lipinski-like', 'Drug-like', or 'Drug-like frag'
        """
        if mode=='default':
            self.setFilterRanges( mode,
                HbDMin=0, HbDMax=99, HbAMin=0, HbAMax=99, MWMin=0, MWMax=9999,
                NatMin=0, NatMax=999, TORSDOFMin=0, TORSDOFMax=32)

        elif mode=='Lipinski-like':
            # http://en.wikipedia.org/wiki/Lipinski%27s_Rule_of_Five
            self.setFilterRanges( mode,
                HbDMin=0, HbDMax=5,HbAMin=0, HbAMax=10, MWMin=0, MWMax=500,
                NatMin=0, NatMax=999, TORSDOFMin=0, TORSDOFMax=32)

        elif mode=='Drug-like':
            # http://en.wikipedia.org/wiki/Lipinski%27s_Rule_of_Five#cite_note-2
            self.setFilterRanges( mode,
                HbDMin=0, HbDMax=5,HbAMin=0, HbAMax=10, MWMin=160, MWMax=480,
                NatMin=20, NatMax=70, TORSDOFMin=0, TORSDOFMax=32)

        elif mode=='Drug-like frag':
            # Values from Fattori's paper
            self.setFilterRanges( mode,
                HbDMin=0, HbDMax=3, HbAMin=0, HbAMax=6, MWMin=160, MWMax=300,
                NatMin=6, NatMax=45, TORSDOFMin=0, TORSDOFMax=32)

        elif mode=='None' or mode==None:
            # Values from Fattori's paper
            self.setFilterRanges( mode,
                HbDMin=None, HbDMax=None, HbAMin=None, HbAMax=None,
                MWMin=None, MWMax=None, NatMin=None, NatMax=None,
                TORSDOFMin=None, TORSDOFMax=None)
        else:
            raise ValueError('BAD FILTER: got %s, expected one of %s'%(
                mode, str(self.filterNames)))


    def filterTable(self, propT, subset=None):
        """
        filter ligands

        LigandList <- FilterLigands.filterTable(propertyTable, subset=None)
        
        is subset is not specified all keys in propertyTables will be
        considered else, only the listed set of keys will be considered for \
        filtering
        """
        from os.path import join
	accepted = []
	rejected = []

        HbAMin = self.HbAMin
        HbAMax = self.HbAMax
        HbDMin = self.HbDMin
        HbDMax = self.HbDMax
        MWMin = self.MWMin
        MWMax = self.MWMax
        NatMin = self.NatMin
        NatMax = self.NatMax
        TORSDOFMin = self.TORSDOFMin
        TORSDOFMax = self.TORSDOFMax

        import os
#        dirname = os.path.basename(os.path.abspath(propT.folder))
        dirname = ""
        table = propT.table
        
        # get all keys in case none were given
        if subset is None:
            subset = table.keys()

        from time import time
        t1 = time()

        if self.filter==None or self.filter=='None':
            accepted = subset
        else:
            # column numbers
            TORS = 0
            HbD  = 1
            HbA  = 2
            MW   = 3
            Nat  = 4
            nstd = 5

            for key in subset:
                v = table[key]
                if self.rejectNStdAt and v[nstd]:
                    rejected.append(key)
                else:
                    lhba =  v[HbA]
                    lhbd =  v[HbD]
                    ltors = v[TORS]
                    lmw =   v[MW]
                    lnat =  v[Nat]
                    if lhba<HbAMin or lhba>HbAMax or \
                       lhbd<HbDMin or lhbd>HbDMax or \
                       ltors<TORSDOFMin or ltors>TORSDOFMax or \
                       lmw<MWMin or lmw>MWMax or \
                       lnat<NatMin or lnat>NatMax:
                        rejected.append(key)
                    else:
                        accepted.append(key)

        self.filteringTime = time()-t1
        at = self.getAtomTypes(propT, accepted)

        # add dirname
        fullkey = [ join(dirname, x) for x in accepted ]
        acceptedList = LigandList( fullkey, at)
        return acceptedList, rejected


    def getAtomTypes(self,  propertyTable, subset=None):

        table = propertyTable.table
        # get all keys in case none were given
        if subset is None:
            subset = table.keys()

        atypes = {}
        
        for key in subset:
            for t in table[key][6]:
               atypes[t] = 1
        return atypes.keys()
    

class LigandList:
    """
    This class represents a list of ligands to be used in a virtual screening
    experiment
    """

    def __init__(self, filename, atomTypes):

        self.filenames = filename # list of filenames containing ligands
        self.atomTypes = atomTypes


    def addAtomTypes(self, atypes):
        for t in atypes:
            if t not in self.atomTypes:
                self.atomTypes.append(t)

if __name__=='__main__':
    import sys, os, time
    print sys.argv[1]
    calcProp = CalculateProperties(verbose=2)
    propT = calcProp.processDirectory(sys.argv[1], "*.pdbqt")

    #propT.saveAsPythonCode('lib1_prop.py')
    propT.saveAsTextFile('lib1.prop')

    # create the PropertyTable from a file
    pt = PropertyTable(filename='lib1.prop')
    
    # create the PropertyTable from a URL
    pt = PropertyTable(url="http://www.scripps.edu/~sanner/collab/lib1.prop")

    lfilter = FilterLigands()
    lfilter.setFilter('None')
    accepted, rejected = lfilter.filterTable(pt)
    print 'No filter in %f second accepted %d rejected %d'%(
        lfilter.filteringTime, len(accepted.filenames), len(rejected))
    print accepted.atomTypes

    lfilter.setFilter('default')
    accepted, rejected = lfilter.filterTable(pt)
    print 'default filter in %f second accepted %d rejected %d'%(
        lfilter.filteringTime, len(accepted.filenames), len(rejected))
    print accepted.atomTypes

    lfilter.setFilter('Lipinski-like')
    accepted, rejected = lfilter.filterTable(pt)
    print 'Lipinski filter in %f second accepted %d rejected %d'%(
        lfilter.filteringTime, len(accepted.filenames), len(rejected))
    print accepted.atomTypes
    
    lfilter.setFilter('Drug-like')
    accepted, rejected = lfilter.filterTable(pt)
    print 'Drug filter in %f second accepted %d rejected %d'%(
        lfilter.filteringTime, len(accepted.filenames), len(rejected))
    print accepted.atomTypes
    
    lfilter.setFilter('Drug-like frag')
    accepted, rejected = lfilter.filterTable(pt)
    print 'Drug frag filter in %f second accepted %d rejected %d'%(
        lfilter.filteringTime, len(accepted.filenames), len(rejected))
    print accepted.atomTypes

##     raise
##     # chaining filters
##     lfilter.setFilter('Lipinski-like')
##     accepted1, rejected1 = lfilter.filterTable(calcProp.propertyTable)

##     lfilter.setFilter('Drug-like frag')
##     accepted, rejected = lfilter.filterTable(calcProp.propertyTable,
##                                              subset=accepted1)
##     print 'Lipinski and Drug frag filter in %f second accepted %d rejected %d'%(
##         lfilter.filteringTime, len(accepted), len(rejected))
##     print lfilter.getAtomTypes(calcProp.propertyTable, accepted)

    
