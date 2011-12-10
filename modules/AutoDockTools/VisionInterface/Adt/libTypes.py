########################################################################
#
# Date: Jan 2006 Authors: Guillaume Vareille, Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
#    Vision Library Types
#
#########################################################################
#
# /home/vareille/.mgltools/1.6.0/Vision/UserLibs/MyDefaultLib/libTypes.py
# Vision will generate this file automatically if it can't find it
#
from LigandDB import LigandDB
from receptor import receptor
from receptor_prepared import receptor_prepared
from dpf_template import dpf_template
from gpf_template import gpf_template
from autogrid_results import autogrid_results
from NetworkEditor.datatypes import AnyArrayType

###################################################
# add new types to your library
###################################################

#class ThingType(AnyArrayType):
#
#    from ThingPackage import Thing
#    def __init__(self, name='thing', color='#995699', shape='rect',
#                 klass=Thing):
#
#        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
#                              klass=klass)
#
## in NetworkEditor.datatypes, you should have a look at the class IntType

class LigandDBType(AnyArrayType):
    
    from AutoDockTools.VisionInterface.Adt.LigandDB import LigandDB
    def __init__(self, name='LigandDB', color='#FFCCFF', shape='rect',
                 klass=LigandDB):

        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
                              klass=klass)

class receptor(AnyArrayType):
    
    from AutoDockTools.VisionInterface.Adt.receptor import receptor
    def __init__(self, name='receptor', color='#99FF33', shape='triangle',
                 klass=receptor):

        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
                              klass=klass)

class receptor_prepared(AnyArrayType):
    
    from AutoDockTools.VisionInterface.Adt.receptor_prepared import receptor_prepared
    def __init__(self, name='receptor_prepared', color='#009900', shape='triangle',
                 klass=receptor_prepared):

        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
                              klass=klass)

class dpf_template(AnyArrayType):
    
    from AutoDockTools.VisionInterface.Adt.dpf_template import dpf_template
    def __init__(self, name='dpf_template', color='#9933FF', shape='triangle',
                 klass=dpf_template):

        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
                              klass=klass)

class gpf_template(AnyArrayType):
    
    from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template
    def __init__(self, name='gpf_template', color='#FF3333', shape='triangle',
                 klass=gpf_template):

        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
                              klass=klass)

class autogrid_results(AnyArrayType):
    
    from AutoDockTools.VisionInterface.Adt.autogrid_results import autogrid_results
    def __init__(self, name='autogrid_results', color='#FF33CC', shape='triangle',
                 klass=autogrid_results):

        AnyArrayType.__init__(self, name=name, color=color, shape=shape, 
                              klass=klass)

## in NetworkEditor.datatypes, you should have a look at the class IntType
