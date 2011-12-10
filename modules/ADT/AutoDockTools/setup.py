#!/usr/bin/env python

from distutils.core import setup
from distutils.command.sdist import sdist
from distutils.command.install_data import install_data
from distutils.command.build_scripts import build_scripts
from glob import glob
import os

########################################################################
# Had to overwrite the prunrefile_list method of sdist to not
# remove automatically the RCS/CVS directory from the distribution.
########################################################################

class modified_sdist(sdist):
    def prune_file_list(self):
        """
        Prune off branches that might slip into the file list as created
        by 'read_template()', but really don't belong there:
          * the build tree (typically 'build')
          * the release tree itself (only an issue if we ran 'sdist
            previously with --keep-temp, or it aborted)
        """
        build = self.get_finalized_command('build')
        base_dir = self.distribution.get_fullname()
        self.filelist.exclude_pattern(None, prefix=build.build_base)
        self.filelist.exclude_pattern(None, prefix=base_dir)

class modified_install_data(install_data):

    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return install_data.run(self)

class modified_build_scripts(build_scripts):
    
    def copy_scripts(self):
        #print "scripts:" , self.scripts
        ## copy runAdt.py to runAdt, inserting "#!python" line at the beginning 
        for script in self.scripts:
            fnew = open(script, "w")
            forig = open(script+".py", "r")
            txt = forig.readlines()
            forig.close()
            fnew.write("#!/usr/bin/env python2.4\n")
            fnew.writelines(txt)
            fnew.close()
        build_scripts.copy_scripts(self)


########################################################################

# list of the python packages to be included in this distribution.
# sdist doesn't go recursively into subpackages so they need to be
# explicitaly listed.
# From these packages only the python modules will be taken
packages = ['AutoDockTools', 
            'AutoDockTools.Tests', 
            'AutoDockTools.Objects',
            #'AutoDockTools.Utilities',
            #'AutoDockTools.Utilities.Tests',
            'AutoDockTools.Utilities24',
            'AutoDockTools.Utilities24.Tests',
            'AutoDockTools.VisionInterface',
            'AutoDockTools.VisionInterface.Adt',
            'AutoDockTools.VisionInterface.Adt.Input',
            'AutoDockTools.VisionInterface.Adt.Macro',
            'AutoDockTools.VisionInterface.Adt.Mapper',
            'AutoDockTools.VisionInterface.Adt.Other',
            'AutoDockTools.VisionInterface.Adt.Output',
            'AutoDockTools.VisionInterface.Adt.Filter',
            'AutoDockTools.CADD',
            'AutoDockTools.CADD.VisionNetworks',
            ]

# list of the python modules not part of a package. Give the path and the
# filename without the extension. i.e you want to add the
# test.py module which is located in MyPack/Tests/ you give
# 'MyPack/Tests/test'
py_modules = ['AutoDockTools/bin/runAdt']

# list of the files that are not python packages but are included in the
# distribution and need to be installed at the proper place  by distutils.
# The list in MANIFEST.in lists is needed for including those files in
# the distribution, data_files in setup.py is needed to install them
# at the right place.
from string import split
data_files = []
def getDataFiles(file_list, directory, names):
    dir_list = split(directory, os.sep)
    if len(dir_list) > 1:
        ff = dir_list[1]
        if ff  == 'Utilities' or ff == 'doc':
            #print "do not include Utilities and doc files"
            return
    
    fs = []
    for name in names:
        ext = os.path.splitext(name)[1]
        #print directory, name, ext, len(ext)
        if ext !=".py" and ext !=".pyc":
            fullname = os.path.join(directory,name)
            if not os.path.isdir(fullname):
                fs.append(fullname)
    if len(fs):
        file_list.append((directory, fs))

os.path.walk("AutoDockTools", getDataFiles, data_files)
#print "data_files:", data_files, len(data_files)


# description of what is going to be included in the distribution and
# installed.
from version import VERSION
setup (name = 'AutoDockTools',
       version = VERSION,
       description = "ADT is a GUI to help set up, launch and analyze AutoDock  dockings",
       author = 'Molecular Graphics Laboratory',
       author_email = 'rhuey@scripps.edu',
       download_url = 'http://www.scripps.edu/~sanner/software/packager.html',
       url = 'http://www.scripps.edu/~sanner/python/adt/autotoolsoverview.html',
       packages = packages,
       py_modules = py_modules,
       data_files = data_files,
       scripts = ['AutoDockTools/bin/runAdt'],
       cmdclass = {'sdist': modified_sdist,
                   'install_data': modified_install_data,
                   'build_scripts': modified_build_scripts
                   },
       )
