## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/APBSParameters.py,v 1.24 2009/09/23 01:34:50 sanner Exp $
#
# $Id: APBSParameters.py,v 1.24 2009/09/23 01:34:50 sanner Exp $
#
import os, types

class Ion:
    """A class that represents ion"""

    def __init__(self, chg = 1.0, con = 1.0, rad = 1.0):
        """Constructor for class Ion"""
        self.charge = chg
        self.concentration = con
        self.radius = rad
        
    def toString(self):
        """Converts to string representation"""
        return '%.3f, %.3f, %.3f' % (self.charge, self.concentration, \
        self.radius)


def checkKeywords(_name, keywords, **kw):
    """test is all kyes in **kw are in list keywords"""
    for key in kw.keys():
        if key not in keywords:
            print 'WARNING: Keyword %s not recognized for %s' % (key, _name)

import numpy
from MolKit.molecule import Atom

class APBSParams:
    """Class to represent APBS parameters.  
       The list of attributes that can be set is found in self.keywords."""

    GAMMA = 0.105 #Apolar energy surface coefficient
    CFAC = 1.5 #Ratio of course grid lengths to molecule bounding box lengths
    MEGABYTES_PER_GRID_POINT = 160.0/1000000.0
    GRID_VALUES = [9, 17, 33, 49, 65, 81, 97, 113, 129, 145, 161,
               177, 193, 209, 225, 241, 257, 273, 289, 305, 321, 337, 353,
            369, 385, 401, 417, 433, 449, 465, 481, 497, 513,
               529, 545, 561, 577, 593, 609, 625, 641, 657, 673, 689]

    FILETYPES = ('', 'OpenDX', 'AVS UCD', 'UHBD')
    CALCULATIONTYPES = ('Solvation energy', 'Binding energy',
                        'Electrostatic potential')
    PBETYPES = ('Linearized', 'Nonlinear')
    BOUNDARYTYPES = ('Zero E', 'Single Debye-Huckel', 'Multiple Debye-Huckel')
    CHARGEDISCRETIZATIONTYPES = ('Trilinear hat-function', 'Cubic B-spline',  'Quintic B-spline')
    SURFACECALCULATIONTYPES = ('No smoothing', 'Harmonic Average',
                               'Cubic B-spline',  '7th Order Polynomial')
    ENERGYOUTPUTTYPES = ('', 'Total', 'Total and Individual')
    FORCEOUTPUTTYPES = ('', 'Total', 'Total and Individual')

    keywords = [
        'APBS_Path',
        'pdb2pqr_Path',
        'pdb2pqr_ForceField', # amber, charm, or parse
        'name',
        'calculationType',
        'pbeType',          
        'boundaryConditions',
        'chargeDiscretization',
        'surfaceCalculation',
        'sdens', #sphere density: see APBS Documentation
        'splineWindow',
        'molecule1Path',
        'molecule2Path',
        'complexPath',
        'energyOutput',
        'forceOutput', 
        'projectFolder',
        'chargeDistributionFile',
        'potentialFile',
        'solventAccessibilityFile',
        'splineBasedAccessibilityFile',
        'VDWAccessibilityFile',
        'ionAccessibilityFile',
        'laplacianOfPotentialFile',
        'energyDensityFile',
        'ionNumberFile',
        'ionChargeDensityFile',
        'xShiftedDielectricFile',
        'yShiftedDielectricFile',
        'zShiftedDielectricFile',
        'kappaFunctionFile',
        'gridPointsX',
        'gridPointsY',
        'gridPointsZ',
        'coarseLengthX',
        'coarseLengthY',
        'coarseLengthZ',
        'coarseCenterX',
        'coarseCenterY',
        'coarseCenterZ',
        'fineLengthX',
        'fineLengthY',
        'fineLengthZ',
        'fineCenterX',
        'fineCenterY',
        'fineCenterZ',
        'proteinDielectric',
        'solventDielectric',
        'solventRadius',  
        'systemTemperature',
        'ions', 
        'saltConcentration',
        'log'
        ]


    def defaultGridSizes(self, mol1, mol2=None):
        """set default fine and corase grid sizes and centers"""
        # get all atomic coordinates
        coords = mol1.findType(Atom).coords
        if self.calculationType == 'Binding energy':
            assert mol2 is not None
            coords += mol2.findType(Atom).coords

        center = (numpy.maximum.reduce(coords)+numpy.minimum.reduce(coords))*0.5
        center = center.tolist()
        self.fineCenterX = self.coarseCenterX = round(center[0],4)
        self.fineCenterY = self.coarseCenterY = round(center[1],4)
        self.fineCenterZ = self.coarseCenterZ = round(center[2],4)

        length = numpy.maximum.reduce(coords) - numpy.minimum.reduce(coords)
        self.coarseLengthX = self.CFAC*(length.tolist())[0] + 10.
        self.coarseLengthY = self.CFAC*(length.tolist())[1] + 10.
        self.coarseLengthZ = self.CFAC*(length.tolist())[2] + 10.
        self.fineLengthX = (length.tolist())[0] + 10.0
        self.fineLengthY = (length.tolist())[1] + 10.0
        self.fineLengthZ = (length.tolist())[2] + 10.0

        
    
    def __init__(self, name='Default', **kw):
        """Constructor for class APBSParams"""
        from mglutil.util.packageFilePath import getBinary, findFilePath, which
        if not hasattr(self,'APBS_Path'):
            try:
                self.APBS_Path = getBinary("apbs", 'binaries')
            except ImportError:
                self.APBS_Path = None
            if not self.APBS_Path:
                self.APBS_Path = which('apbs')
                 
        self.pdb2pqr_Path = findFilePath('pdb2pqr.py','MolKit.pdb2pqr')
        self.name = name    # name of the APBS Param instance
        
        self.pdb2pqr_ForceField = 'amber'
        self.calculationType = 'Electrostatic potential'
        self.pbeType = 'Linearized'
        self.boundaryConditions = 'Single Debye-Huckel'
        self.chargeDiscretization = 'Cubic B-spline'
        self.surfaceCalculation = 'Cubic B-spline'
        self.sdens = 10.0
        self.splineWindow = 0.3
        self.molecule1Path = ''
        self.molecule2Path = ''
        self.complexPath = ''
        self.energyOutput = 'Total'
        self.forceOutput = ''
        self.projectFolder = 'apbs-project'
        self.chargeDistributionFile = ''
        self.potentialFile = 'OpenDX'
        self.solventAccessibilityFile = ''
        self.splineBasedAccessibilityFile = ''
        self.VDWAccessibilityFile = ''
        self.ionAccessibilityFile = ''
        self.laplacianOfPotentialFile = ''
        self.energyDensityFile = ''
        self.ionNumberFile = ''
        self.ionChargeDensityFile = ''
        self.xShiftedDielectricFile = ''
        self.yShiftedDielectricFile = ''
        self.zShiftedDielectricFile = ''
        self.kappaFunctionFile = ''

        #Grid parameters
        self.gridPointsX = 65
        self.gridPointsY = 65
        self.gridPointsZ = 65
        self.coarseLengthX = 40
        self.coarseLengthY = 50
        self.coarseLengthZ = 60
        self.coarseCenterX = 0
        self.coarseCenterY = 0
        self.coarseCenterZ = 0
        self.fineLengthX = 20
        self.fineLengthY = 35
        self.fineLengthZ = 30
        self.fineCenterX = 0
        self.fineCenterY = 0
        self.fineCenterZ = 0

        # Physical Parameters
        self.proteinDielectric = 2.0
        self.solventDielectric = 78.54
        self.solventRadius = 1.4
        self.systemTemperature = 298.15
        self.saltConcentration = 0.01
        self.ions = []

        apply( self.Set, (), kw )

    def Set(self, check=1, **kw):
        """Sets APBSParams member variable(s)"""
        if check:
            apply( checkKeywords, ('APBSParam object'+self.name,
                                   self.keywords), kw)

        val = kw.get('APBS_Path', None)
        if val:
            if os.access(val, os.X_OK):
                self.APBS_Path = val
        # either we have pdb2pqr under MolKit and independent of APBS
        # then we can simply run the python code in this interpreter
        # else we HAVE TO set the path to a place where pdb2qpr is
        # else we refues non pqr files as valid input
        val = kw.get('pdb2pqr_Path', None)
        if val:
            if os.path.exists(val):
                self.pdb2pqr_Path = val
            
        val = kw.get('name', None)
        if val:
            self.name = val

        val = kw.get("pdb2pqr_ForceField", None)
        if val:
            assert val in ['amber', 'charmm', 'parse']
            self.pdb2pqr_ForceField = val
        
        val = kw.get("projectFolder", None)
        if val:
            self.projectFolder = val

        val = kw.get('molecule1Path', None)
        if val:
            self.molecule1Path = val

        val = kw.get('molecule2Path', None)
        if val:
            self.molecule2Path = val

        val = kw.get('complexPath', None)
        if val:
            self.complexPath = val
        
        #Calculation parameters 

        val = kw.get('calculationType', None)
        if val:
            assert val in self.CALCULATIONTYPES
            if val == 'Solvation energy' or val == 'Electrostatic potential':
                assert self.molecule1Path is not None
                self.calculationType = val
            elif val == 'Binding energy':
                assert self.molecule1Path is not None
                assert self.molecule2Path is not None
                assert self.complexPath is not None
            self.calculationType = val
            
        val = kw.get('pbeType', None)            
        if val:
            assert val in self.PBETYPES
            self.pbeType = val

        val = kw.get('boundaryConditions', None)            
        if val:
            assert val in self.BOUNDARYTYPES
            self.boundaryConditions = val
            
        val = kw.get('chargeDiscretization', None)            
        if val:
            #assert val in self.CHARGEDISCRETIZATIONTYPES
            self.chargeDiscretization = val

        val = kw.get('surfaceCalculation', None)            
        if val:
            #assert val in self.SURFACECALCULATIONTYPES
            self.surfaceCalculation = val

        val = kw.get('sdens', None)            
        if val:
            assert val != '' and (isinstance(val, types.FloatType) or \
                                isinstance(val, types.IntType)) and val > 0.0
            self.splineWindow = val

        val = kw.get('splineWindow', None)            
        if val:
            assert val != '' and (isinstance(val, types.FloatType) or \
                                isinstance(val, types.IntType)) and val > 0.0
            self.splineWindow = val
            
        val = kw.get('energyOutput', None)            
        if val:
            assert val in self.ENERGYOUTPUTTYPES
            self.energyOutput = val

        val = kw.get('forceOutput', None)            
        if val:
            assert val in self.FORCEOUTPUTTYPES
            self.forceOutput = val

        val = kw.get('chargeDistributionFile', None)            
        if val:
            assert val in self.FILETYPES
            self.chargeDistributionFile = val
            
        val = kw.get('chargeDistributionFile', None)            
        if val:
            assert val in self.FILETYPES
            self.chargeDistributionFile = val

        val = kw.get('potentialFile', None)            
        if val:
            assert val in self.FILETYPES
            self.potentialFile = val

        val = kw.get('solventAccessibilityFile', None)            
        if val:
            assert val in self.FILETYPES
            self.solventAccessibilityFile = val
            
        val = kw.get('splineBasedAccessibilityFile', None)            
        if val:
            assert val in self.FILETYPES
            self.splineBasedAccessibilityFile = val

        val = kw.get('VDWAccessibilityFile', None)            
        if val:
            assert val in self.FILETYPES
            self.VDWAccessibilityFile = val

        val = kw.get('ionAccessibilityFile', None)            
        if val:
            assert val in self.FILETYPES
            self.ionAccessibilityFile = val

        val = kw.get('laplacianOfPotentialFile', None)            
        if val:
            assert val in self.FILETYPES
            self.laplacianOfPotentialFile = val

        val = kw.get('energyDensityFile', None)            
        if val:
            assert val in self.FILETYPES
            self.energyDensityFile = val

        val = kw.get('ionNumberFile', None)            
        if val:
            assert val in self.FILETYPES
            self.ionNumberFile = val

        val = kw.get('ionChargeDensityFile', None)            
        if val:
            assert val in self.FILETYPES
            self.ionChargeDensityFile = val

        val = kw.get('xShiftedDielectricFile', None)            
        if val:
            assert val in self.FILETYPES
            self.xShiftedDielectricFile = val

        val = kw.get('yShiftedDielectricFile', None)            
        if val:
            assert val in self.FILETYPES
            self.yShiftedDielectricFile = val

        val = kw.get('zShiftedDielectricFile', None)            
        if val:
            assert val in self.FILETYPES
            self.zShiftedDielectricFile = val

        val = kw.get('kappaFunctionFile', None)            
        if val:
            assert val in self.FILETYPES
            self.kappaFunctionFile = val

        #Grid parameters
        val = kw.get('gridPointsX', None)            
        if val:
            assert val in self.GRID_VALUES
            self.gridPointsX = val

        val = kw.get('gridPointsY', None)            
        if val:
            assert val in self.GRID_VALUES
            self.gridPointsY = val

        val = kw.get('gridPointsZ', None)            
        if val:
            assert val in self.GRID_VALUES
            self.gridPointsZ = val

        val = kw.get('coarseLengthX', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.coarseLengthX = val

        val = kw.get('coarseLengthY', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.coarseLengthY = val

        val = kw.get('coarseLengthZ', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.coarseLengthZ = val

        val = kw.get('coarseCenterX', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.coarseCenterX = val

        val = kw.get('coarseCenterY', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.coarseCenterY = val

        val = kw.get('coarseCenterZ', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.coarseCenterZ = val

        val = kw.get('coarseCenterZ', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.coarseCenterZ = val

        val = kw.get('fineLengthX', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.fineLengthX = val

        val = kw.get('fineLengthY', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.fineLengthY = val

        val = kw.get('fineLengthZ', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.fineLengthZ = val

        val = kw.get('fineCenterX', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.fineCenterX = val

        val = kw.get('fineCenterY', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.fineCenterY = val

        val = kw.get('fineCenterZ', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.fineCenterZ = val

        #Physical parameters       
        val = kw.get('proteinDielectric', None)            
        if val:
            assert isinstance(val, types.FloatType)
            self.proteinDielectric = val

        val = kw.get('solventDielectric', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))
            self.solventDielectric = val

        val = kw.get('solventRadius', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.solventRadius = val

        val = kw.get('systemTemperature', None)            
        if val:
            assert (isinstance(val, types.FloatType) or isinstance(val, types.IntType))\
                                                                   and val > 0.0
            self.systemTemperature = val

        val = kw.get('ions', [])            
        for ion in val:
            assert ion.__doc__ == Ion.__doc__
            self.ions.append(ion)


    def apbsWriteCalculationParams(self, fp, molname):
        """None <--- apbsWriteCalculationParams(fp)\n
           Writes APBS Calculation Parameters into fp\n
        """
        if(self.pbeType=='Linearized'):
            fp.write('\tlpbe\n')
        else:
            fp.write('\tnpbe\n')
            
        if(self.boundaryConditions=='Zero E'):
            fp.write('\tbcfl zero\n')
        elif(self.boundaryConditions=='Single Debye-Huckel'):
            fp.write('\tbcfl sdh\n')
        else: fp.write('\tbcfl mdh\n')

        if(self.chargeDiscretization=='Trilinear hat-function'):
            fp.write('\tchgm spl0\n')
        elif self.chargeDiscretization == 'Cubic B-spline':
            fp.write('\tchgm spl2\n')
        else:
            fp.write('\tchgm spl4\n')
        
        
        if(self.surfaceCalculation=='No smoothing'):
            fp.write('\tsrfm mol\n')
            fp.write('\tsdens %.3f\n'%(self.sdens))
        elif(self.surfaceCalculation=='Harmonic Average'):
            fp.write('\tsrfm smol\n')
            fp.write('\tsdens %.3f\n'%(self.sdens))
        elif self.surfaceCalculation == 'Cubic B-spline':
            fp.write('\tsrfm spl2\n')
            fp.write('\tswin %.3f\n'%(self.splineWindow))
        else:
            fp.write('\tsrfm spl4\n')
            fp.write('\tswin %.3f\n'%(self.splineWindow))
        
                        
        if(self.energyOutput==''):
            fp.write('\tcalcenergy no\n')
        elif(self.energyOutput=='Total'):
            fp.write('\tcalcenergy total\n')
        else: fp.write('\tcalcenergy comps\n')

        if(self.forceOutput==''):
            fp.write('\tcalcforce no\n')
        elif(self.forceOutput=='Total'):
            fp.write('\tcalcforce total\n')
        else: fp.write('\tcalcforce comps\n')

        tempFileString = molname + '.chargeDistribution'
        if  (self.chargeDistributionFile=='OpenDX'):    
            fp.write('\twrite charge dx %s\n'  % tempFileString)
        elif(self.chargeDistributionFile=='AVS UCD'):
            fp.write('\twrite charge avs %s\n' % tempFileString)
        elif(self.chargeDistributionFile=='UHBD'):
            fp.write('\twrite charge uhbd %s\n'%tempFileString)

        tempFileString = molname +'.potential'
        if  (self.potentialFile=='OpenDX'):
            fp.write('\twrite pot dx %s\n'  % tempFileString)
        elif(self.potentialFile=='AVS UCD'):
            fp.write('\twrite pot avs %s\n' % tempFileString)
        elif(self.potentialFile=='UHBD'):
            fp.write('\twrite pot uhbd %s\n'%tempFileString)

        tempFileString = molname + '.solventAccessibility'
        if  (self.solventAccessibilityFile=='OpenDX'):
            fp.write('\twrite smol dx %s\n'  % tempFileString)
        elif(self.solventAccessibilityFile=='AVS UCD'):
            fp.write('\twrite smol avs %s\n' % tempFileString)
        elif(self.solventAccessibilityFile=='UHBD'):
            fp.write('\twrite smol uhbd %s\n'%tempFileString)

        tempFileString = molname + '.splineBasedAccessibility'
        if (self.splineBasedAccessibilityFile=='OpenDX'):
            fp.write('\twrite sspl dx %s\n'  % tempFileString)
        elif(self.splineBasedAccessibilityFile=='AVS UCD'):
            fp.write('\twrite sspl avs %s\n' % tempFileString)
        elif(self.splineBasedAccessibilityFile=='UHBD'):
            fp.write('\twrite sspl uhbd %s\n'% tempFileString)

        tempFileString = molname + '.VDWAccessibility'
        if  (self.VDWAccessibilityFile=='OpenDX'): 
            fp.write('\twrite vdw dx %s\n'  % tempFileString)
        elif(self.VDWAccessibilityFile=='AVS UCD'):
            fp.write('\twrite vdw avs %s\n' % tempFileString)
        elif(self.VDWAccessibilityFile=='UHBD'):
            fp.write('\twrite vdw uhbd %s\n'% tempFileString)

        tempFileString = molname + '.ionAccessibility'
        if  (self.ionAccessibilityFile=='OpenDX'):
            fp.write('\twrite ivdw dx %s\n'  % tempFileString)
        elif(self.ionAccessibilityFile=='AVS UCD'):
            fp.write('\twrite ivdw avs %s\n' % tempFileString)
        elif(self.ionAccessibilityFile=='UHBD'):
            fp.write('\twrite ivdw uhbd %s\n'% tempFileString)

        tempFileString = molname + '.laplacianOfPotential'
        if  (self.laplacianOfPotentialFile=='OpenDX'):
            fp.write('\twrite lap dx %s\n'  % tempFileString)
        elif(self.laplacianOfPotentialFile=='AVS UCD'):
            fp.write('\twrite lap avs %s\n' % tempFileString)
        elif(self.laplacianOfPotentialFile=='UHBD'):
            fp.write('\twrite lap uhbd %s\n'% tempFileString)

        tempFileString = molname + '.energyDensity'
        if  (self.energyDensityFile=='OpenDX'): 
            fp.write('\twrite edens dx %s\n'  % tempFileString)
        elif(self.energyDensityFile=='AVS UCD'): 
            fp.write('\twrite edens avs %s\n' % tempFileString)
        elif(self.energyDensityFile=='UHBD'):
            fp.write('\twrite edens uhbd %s\n'% tempFileString)

        tempFileString = molname +'.ionNumber'
        if  (self.ionNumberFile=='OpenDX'):
            fp.write('\twrite ndens dx %s\n'  % tempFileString)
        elif(self.ionNumberFile=='AVS UCD'): 
            fp.write('\twrite ndens avs %s\n' % tempFileString)
        elif(self.ionNumberFile=='UHBD'): 
            fp.write('\twrite ndens uhbd %s\n'% tempFileString)

        tempFileString = molname + '.ionChargeDensity'
        if  (self.ionChargeDensityFile=='OpenDX'):
            fp.write('\twrite qdens dx %s\n'  % tempFileString)
        elif(self.ionChargeDensityFile=='AVS UCD'):
            fp.write('\twrite qdens avs %s\n' % tempFileString)
        elif(self.ionChargeDensityFile=='UHBD'):
            fp.write('\twrite qdens uhbd %s\n'% tempFileString)

        tempFileString = molname + '.xShiftedDielectric'
        if  (self.xShiftedDielectricFile=='OpenDX'):
            fp.write('\twrite dielx dx %s\n'  % tempFileString)
        elif(self.xShiftedDielectricFile=='AVS UCD'):
            fp.write('\twrite dielx avs %s\n' % tempFileString)
        elif(self.xShiftedDielectricFile=='UHBD'):
            fp.write('\twrite dielx uhbd %s\n'% tempFileString)

        tempFileString = molname + '.yShiftedDielectric'
        if  (self.yShiftedDielectricFile=='OpenDX'):
            fp.write('\twrite diely dx %s\n'  % tempFileString)
        elif(self.yShiftedDielectricFile=='AVS UCD'):
            fp.write('\twrite diely avs %s\n' % tempFileString)
        elif(self.yShiftedDielectricFile=='UHBD'):
            fp.write('\twrite diely uhbd %s\n'%tempFileString)

        tempFileString = molname + '.zShiftedDielectric'
        if  (self.zShiftedDielectricFile=='OpenDX'): 
            fp.write('\twrite dielz dx %s\n'  % tempFileString)
        elif(self.zShiftedDielectricFile=='AVS UCD'):
            fp.write('\twrite dielz avs %s\n' % tempFileString)
        elif(self.zShiftedDielectricFile=='UHBD'):
            fp.write('\twrite dielz uhbd %s\n'%tempFileString)

        tempFileString = molname + '.kappaFunction'
        if  (self.kappaFunctionFile=='OpenDX'):
            fp.write('\twrite kappa dx %s\n'  % tempFileString)
        elif(self.kappaFunctionFile=='AVS UCD'):
            fp.write('\twrite kappa avs %s\n' % tempFileString)
        elif(self.kappaFunctionFile=='UHBD'):
            fp.write('\twrite kappa uhbd %s\n'%tempFileString)
        fp.write('\n')

    def apbsWriteGridParams(self, fp):
        """None <--- apbsWriteGridParams(fp)\n
           Writes APBS Grid Parameters into fp\n
        """
        fp.write('\tdime %d %d %d\n\n'%(
            self.gridPointsX,self.gridPointsY, self.gridPointsZ))
        fp.write('\tcglen %.3f %.3f %.3f\n'%(
            self.coarseLengthX,self.coarseLengthY, self.coarseLengthZ))
        fp.write('\tcgcent %.3f %.3f %.3f\n'%(
            self.coarseCenterX,self.coarseCenterY, self.coarseCenterZ))
        fp.write('\tfglen %.3f %.3f %.3f\n'%(
            self.fineLengthX,self.fineLengthY, self.fineLengthZ))
        fp.write('\tfgcent %.3f %.3f %.3f\n'%(
            self.fineCenterX,self.fineCenterY, self.fineCenterZ))
        fp.write('\n')


    def apbsWritePhysicsParams(self, fp):
        """None <--- apbsWritePhysicsParams(fp)\n
           Writes APBS Physics Parameters into fp\n
        """
        #fp.write('\tgamma %.3f\n'%(self.GAMMA)) # NOTE: CONSTANT
        fp.write('\ttemp %.3f\n'%(self.systemTemperature))
        fp.write('\tsrad %.3f\n'%(self.solventRadius))
        fp.write('\tsdie %.3f\n'%(self.solventDielectric))
        fp.write('\tpdie %.3f\n'%(self.proteinDielectric))
        for i in range(0, len(self.ions)):
            fp.write('\tion %s\n'%(self.ions[i].toString()))
        if self.saltConcentration:
            fp.write('\tion 1.000, %.3f, 2.000\n'%(self.saltConcentration))
            fp.write('\tion -1.000, %.3f, 2.000\n'%(self.saltConcentration))
        fp.write('\n')


    def apbsWriteSolvationEnergy(self, fp):
        """None <--- apbsWriteSolvationEnergy(fp)\n
           Writes APBS Solvation Energy Parameters into fp\n
        """
        fp.write('READ\n')
        fp.write('\tmol pqr %s\n'%(self.molecule1Path))
        fp.write('END\n\n')
        fp.write('ELEC\n')
        fp.write('\tmg-auto\n')
        fp.write('\tmol 1\n')
        file_name, ext = os.path.splitext(self.molecule1Path)
        mol_name = os.path.split(file_name)[-1]
        self.apbsWriteCalculationParams(fp, mol_name)
        self.apbsWriteGridParams(fp)
        self.apbsWritePhysicsParams(fp)
        fp.write('END\n\n')
        fp.write('ELEC\n')
        fp.write('\tmg-auto\n')
        fp.write('\tmol 1\n')
        self.apbsWriteCalculationParams(fp, mol_name + '_Vacuum')
        self.apbsWriteGridParams(fp)
        tempSolventDielectric = self.solventDielectric
        self.solventDielectric = 1.0
        tempIons = self.ions
        tempSaltConcentration = self.saltConcentration
        self.ions = []
        self.saltConcentration = None
        self.apbsWritePhysicsParams(fp)
        self.solventDielectric = tempSolventDielectric
        self.ions = tempIons
        self.saltConcentration = tempSaltConcentration
        fp.write('END\n\n')
        fp.write('PRINT\n')
        fp.write('\telecEnergy 1 - 2\n')
        fp.write('END\n\n')
        fp.write('QUIT\n')

    def apbsWriteBindingEnergy(self, fp):
        """None <--- apbsWriteBindingEnergy(fp)\n
           Writes APBS Binding Energy Parameters into fp\n
        """
        fp.write('READ\n')
        fp.write('\tmol pqr %s\n'%(self.molecule1Path))
        fp.write('\tmol pqr %s\n'%(self.molecule2Path))
        fp.write('\tmol pqr %s\n'%(self.complexPath))
        fp.write('END\n\n')
        fp.write('ELEC\n')
        fp.write('\tmg-auto\n')
        fp.write('\tmol 1\n')
        file_name, ext = os.path.splitext(self.molecule1Path)
        mol_name = os.path.split(file_name)[-1]        
        self.apbsWriteCalculationParams(fp, mol_name)
        self.apbsWriteGridParams(fp)
        self.apbsWritePhysicsParams(fp)
        fp.write('END\n\n')
        fp.write('ELEC\n')
        fp.write('\tmg-auto\n')
        fp.write('\tmol 2\n')
        file_name, ext = os.path.splitext(self.molecule2Path)
        mol_name = os.path.split(file_name)[-1]
        self.apbsWriteCalculationParams(fp, mol_name)
        self.apbsWriteGridParams(fp)
        self.apbsWritePhysicsParams(fp)
        fp.write('END\n\n')
        fp.write('ELEC\n')
        fp.write('\tmg-auto\n')
        fp.write('\tmol 3\n')
        file_name, ext = os.path.splitext(self.complexPath)
        mol_name = os.path.split(file_name)[-1]
        self.apbsWriteCalculationParams(fp, mol_name)
        self.apbsWriteGridParams(fp)
        self.apbsWritePhysicsParams(fp)
        fp.write('END\n\n')
        fp.write('PRINT\n')
        fp.write('\telecEnergy 3 - 2 - 1\n')
        fp.write('END\n\n')
        fp.write('QUIT\n')

    def apbsWriteElectrostaticPotential(self, fp):
        """None <--- apbsWriteElectrostaticPotential(fp)\n
           Writes APBS Electrostatic Potential Parameters into fp\n
        """
        fp.write('READ\n')
        fp.write('\tmol pqr %s\n'%(self.molecule1Path))
        fp.write('END\n\n')
        fp.write('ELEC\n')
        fp.write('\tmg-auto\n')
        fp.write('\tmol 1\n')
        file_name, ext = os.path.splitext(self.molecule1Path)
        mol_name = os.path.split(file_name)[-1]        
        self.apbsWriteCalculationParams(fp, mol_name)
        self.apbsWriteGridParams(fp)
        self.apbsWritePhysicsParams(fp)
        fp.write('END\n\n')
        fp.write('PRINT\n')
        fp.write('\telecEnergy 1\n')
        fp.write('END\n\n')
        fp.write('QUIT\n')

    def SaveAPBSInput(self, filename):
        """None <--- apbsWriteElectrostaticPotential(filename)\n
           Saves APBS Input Parameters in a file named filename \n
        """
        fp = open(filename, 'wb+')
        if(self.calculationType=='Solvation energy'):
            self.apbsWriteSolvationEnergy(fp)
        elif(self.calculationType=='Binding energy'):
            self.apbsWriteBindingEnergy(fp)
        else: self.apbsWriteElectrostaticPotential(fp)
        fp.close()
