########################################################################
#
# Copyright: Sargis Dallakyan (sargis@scripps.edu)
#            Michel Sanner (sanner@scripps.edu)
#            The Scripps Research Institute (TSRI)
#            Molecular Graphics Lab
#            La Jolla, CA 92037, USA
#
#########################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/preferences.py,v 1.10 2010/11/11 21:48:34 sargis Exp $
#
# $Id: preferences.py,v 1.10 2010/11/11 21:48:34 sargis Exp $
#
"""
This module contains UserPreference class that is used to store user preferences and settings.
"""

import os
import pickle
import types
from UserDict import UserDict
from mglutil.util.packageFilePath import getResourceFolderWithVersion

class UserPreference(UserDict):
    """
    Class to let the user define Preferences.
    a preference is made of a name, a current value, a possibly empty list
    of valid values.
    preferences can be added using the add method
    and set using the set method
    """

    def __init__(self, ):
        UserDict.__init__(self)
        self.dirty = 0 # used to remember that something changed
        self.resourceFile = None
        resourceFolder = getResourceFolderWithVersion()
        if resourceFolder is None:
            return
        self.resourceFile =  os.path.join(resourceFolder, '.settings')
        self.defaults = {}
        self.settings = {}
        if os.path.exists(self.resourceFile):
            try:
                pkl_file = open(self.resourceFile)
                self.settings = pickle.load(pkl_file)
                pkl_file.close()
            except Exception, inst:
                print inst, "Error in ", __file__
                
    def add(self, name, value, validValues = [], validateFunc=None,
            callbackFunc=[], doc='', category="General"):
        """add a userpreference. A name and a value are required,
        a list of valide values can be provoded as well as a function that
        can be called to validate the value. A callback function can be
        specified, it will be called when the value is set with the old value
        and the new value passed as an argument"""

#        if name in self.data.keys():
#            # doesn't create the userpreference if the name already exists:
#            return
        if len(validValues):
            assert value in validValues
        if validateFunc:
            assert callable(validateFunc)
        if callbackFunc != []:
            assert type(callbackFunc) is types.ListType and \
                   len(filter(lambda x: not callable(x), callbackFunc))==0
        

        self[name] = { 'value':value, 'validValues':validValues,
                       'validateFunc': validateFunc,
                       'callbackFunc': callbackFunc,
                       'doc':doc ,
                       'category':category}
        self.set(name, value)
        self.dirty = 1
        

    def set(self, name, value):
        if not self.data.has_key(name):
            self.settings[name] = value
            return
        if self.resourceFile is None:
            return
        self.settings[name] = value
        entry = self.data[name]
        try:
            if entry.has_key('validValues') and len(entry['validValues']):
                if not value in entry['validValues']:
                    #msg = " is not a valid value, value has to be in %s" % str(entry['validValues'])
                    #print value, msg
                    return
            if entry.has_key('validateFunc') and entry['validateFunc']:
                if not entry['validateFunc'](value):
                    msg = " is not a valid value, try the Info button"
                    #print value, msg
                    return
        except Exception, inst:
            print __file__, inst

        oldValue = entry['value']
        entry['value'] = value
        if entry['callbackFunc']!=[]:
            for cb in  entry['callbackFunc']:
                cb(name,oldValue, value)
        self.dirty = 1
 
    def addCallback(self, name, func):
        assert callable(func)
        assert self.data.has_key(name)
        entry = self.data[name]
        entry['callbackFunc'].append(func)

    def removeCallback(self, name, func):
        assert self.data.has_key(name) and \
               func in self.data[name]['callbackFunc']
        entry = self.data[name]
        entry['callbackFunc'].remove(func)
        
        
    def save(self, filename):
        """save the preferences to a file"""
        pass
        self.dirty = 0 # clean now !


    def loadSettings(self):
        if self.resourceFile is None:
            return
        settings = {}        
        if os.path.exists(self.resourceFile):
            pkl_file = open(self.resourceFile)
            settings = pickle.load(pkl_file)
            pkl_file.close()
        for key, value in settings.items():
            self.set(key, value)
        
    def saveAllSettings(self):
        output = open(self.resourceFile, 'w')
        pickle.dump(self.settings, output)
        output.close()
    
    def saveSingleSetting(self, name, value):
        if os.path.exists(self.resourceFile):
            pkl_file = open(self.resourceFile)
            settings = pickle.load(pkl_file)
            pkl_file.close()
        else:
            settings = {}
        settings[name] = value
        output = open(self.resourceFile, 'w')
        pickle.dump(settings, output)
        output.close()
        
    def saveDefaults(self):
        if self.resourceFile is None:
            return
        for key, value in self.data.items():
            self.defaults[key] = value
    
    def restoreDefaults(self):
        for key, value in self.defaults.items():
            self.set(key, value)
    