'''
This script is described at: http://www.pymolwiki.org/index.php/plugindirectory

PyMOL personal plugin directory
Recommended location of this file: ~/xxx/Pymol-script-repo/plugins/__init__.py

### Make these settings in ~/.pymolrc ####################
import sys
sys.path.append('C:/Users/YOURNAME/Documents/Pymol-script-repo')
import plugins
'''
 
import os, sys, traceback
 
# import pymolplugins (allow different name)
pymolplugins = sys.modules[__name__]
 
import pmg_tk.PMGApp
x__initializePlugins = pmg_tk.PMGApp.initializePlugins
 
def initializePlugins(self):
    '''
    Overloaded version of pmg_tk.PMGApp.initializePlugins
    See pmg_tk/PMGApp.py
    '''
    # load global plugins
    x__initializePlugins(self)
 
    # load user plugins
    modules = set()
    for path in pymolplugins.__path__:
        for filename in os.listdir(path):
            name, _, ext = filename.partition('.')
            if ext not in ['py', 'pyc', 'pyo']:
                if os.path.isdir(os.path.join(path, filename)):
                    modules.add(filename)
            elif name != '__init__':
                modules.add(name)
    for name in modules:
        mod_name = pymolplugins.__name__ + '.' + name
        try:
            __import__(mod_name, level=0)
            mod = sys.modules[mod_name]
            if hasattr(mod,'__init_plugin__'):
                mod.__init_plugin__(self)
            elif hasattr(mod,'__init__'):
                mod.__init__(self)
        except:
            print "Exception in plugin '%s' -- Traceback follows..."%name
            traceback.print_exc()
            print "Error: unable to initialize plugin '%s'."%name
 
# overload method
pmg_tk.PMGApp.initializePlugins = initializePlugins
