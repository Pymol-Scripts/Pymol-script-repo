"""Impements Recent Files menues"""
# $Header: /opt/cvs/python/packages/share1.5/mglutil/util/recentFiles.py,v 1.16 2010/10/05 17:25:40 sargis Exp $
#
# $Id: recentFiles.py,v 1.16 2010/10/05 17:25:40 sargis Exp $
import os, pickle
from mglutil.util.packageFilePath import getResourceFolderWithVersion
import Tkinter

class RecentFiles:
    """Class to store Recent Files"""
    def __init__(self, masterApp, masterMenu, filePath=None, 
                 menuLabel="Open recent", underline=5, index=0):
        """Construct recent files categories. If filePath is not provided
        mglutil/recent.pkl is used to store and load the data"""
        if not filePath: #use "mglutil/recent.pkl" to store the data
            filePath = getResourceFolderWithVersion()
            if filePath is None:
                return
            filePath += os.sep + "mglutil" + os.sep + "recent.pkl"
        if os.path.exists(filePath):
            try:
                self.categories = pickle.load(open(filePath))
            except Exception, inst:
                #print inst
                #print "Couldn't Load Recent Files."
                self.categories = {}
        else:
            self.categories = {}
        self.resourceFilePath = filePath
        self.checkCategories()
        if masterMenu != None:
            self.gui(masterMenu, menuLabel, underline=underline, index=index)
        else:
            self.mainMenu = None
        self.masterApp = masterApp
        
    def checkCategories(self):
        """Loops through self.categories to check if recent file still exists.
        If not, removes the file from the list"""
        for category in self.categories.keys():
            newList = [x for x in self.categories[category] if os.path.exists(x[0])]
            if len(newList):
                self.categories[category] = newList
            else:
                self.categories.pop(category)
                                
    def add(self, filePath, cmdStr, category="Documents"):
        """Add file to self.categories[category] list. 
        First element in this list is the file. 
        Second is the command string - cmdStr."""
        if not filePath:
            return
        if hasattr(self, 'categories') is False: 
            return

        if not self.categories.has_key(category):
            self.categories[category] = []
            #self.menuList[category] = Tkinter.Menu(self.mainMenu)
            #self.mainMenu.add_cascade(label=category, 
            #                          menu=self.menuList[category])
        if os.path.exists(filePath):
            filePath = os.path.abspath(filePath)
            if [filePath,cmdStr] in self.categories[category]:
                index = self.categories[category].index([filePath,cmdStr])
                self.categories[category].pop(index)
                if self.mainMenu != None : self.mainMenu.delete(index+1)
            if len(self.categories[category]) > 20:
                self.categories[category].pop()
                #self.menuList[category].delete(10,Tkinter.END)  
            self.categories[category].insert(0, [filePath,cmdStr])
            if self.mainMenu != None : self.mainMenu.insert(0,'command', label=filePath, 
                                        command=self.callback([filePath,cmdStr]))
        self.dumpCategories()
        
    def dumpCategories(self, filePath=None):
        """Calls pickle.dump(self.categories, open(filePath,'w'))"""
        if not filePath:
            filePath = self.resourceFilePath
        try:
            pickle.dump(self.categories, open(filePath,'w'))
        except Exception, inst:
            print "Failed to save recent files"
            print inst
        
    def gui(self, masterMenu, menuLabel, underline=None, index=0):
        self.mainMenu = Tkinter.Menu(masterMenu)
        masterMenu.insert_cascade(index, label=menuLabel,
                                  menu=self.mainMenu, underline=underline)
        self.menuList = {}
        for category in self.categories:
#             self.menuList[category] = Tkinter.Menu(self.mainMenu)
#             self.mainMenu.add_cascade(label=category, 
#                                       menu=self.menuList[category])
             for listItem in self.categories[category]:
                 self.mainMenu.add_command(label=listItem[0], 
                                           command=self.callback(listItem))

    def callback(self, listItem):
        def call(listItem=listItem):
            masterApp = self.masterApp
            eval("masterApp." + listItem[1] + '(' + repr(listItem[0]) +')') 
        return call
