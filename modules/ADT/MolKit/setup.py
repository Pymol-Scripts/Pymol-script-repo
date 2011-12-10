#!/usr/bin/env python

from distutils.core import setup
from distutils.command.sdist import sdist
from distutils.command.install_data import install_data
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
########################################################################

# list of the python packages to be included in this distribution.
# sdist doesn't go recursively into subpackages so they need to be
# explicitaly listed.
# From these packages only the python modules will be taken
packages = ['MolKit', 'MolKit.data',
            'MolKit.Tests', 'MolKit.VisionInterface',
            'MolKit.VisionInterface.Tests',
            'MolKit.VisionInterface.Tests.Data',
            'MolKit.pdb2pqr', 'MolKit.pdb2pqr.src',
            'MolKit.pdb2pqr.extensions', 'MolKit.pdb2pqr.pdb2pka',
            'MolKit.pdb2pqr.propka', 'MolKit.pdb2pqr.pdb2pka.ligandclean',
            'MolKit.pdb2pqr.pdb2pka.substruct']

# list of the python modules not part of a package. Give the path and the
# filename without the extension. i.e you want to add the
# test.py module which is located in MyPack/Tests/ you give
# 'MyPack/Tests/test'
py_modules = []

# list of the files that are not python packages but are included in the
# distribution and need to be installed at the proper place  by distutils.
# The list in MANIFEST.in lists is needed for including those files in
# the distribution, data_files in setup.py is needed to install them
# at the right place.
data_files = []
## for dir in ['MolKit/data','MolKit',
##             'MolKit/Tests/Data','MolKit/VisionInterface/Tests',
##             'MolKit/CVS', 'MolKit/data/CVS',
##             'MolKit/Tests/CVS', 'MolKit/VisionInterface/CVS',
##             'MolKit/VisionInterface/Tests/CVS','MolKit/VisionInterface/Tests'
##             'MolKit/Tests/Data/CVS',
##             'MolKit/pdb2pqr','MolKit/pdb2pqr/doc']:
##     files = []
##     for f in glob(os.path.join(dir, '*')):
##         if f[-3:] != '.py' and f[-4:-1] != '.py' and os.path.isfile(f):
##             files.append(f)
##     data_files.append((dir, files))

def getDataFiles(file_list, directory, names):
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

os.path.walk("MolKit", getDataFiles, data_files)


# description of what is going to be included in the distribution and
# installed.
from version import VERSION
setup (name = 'MolKit',
       version = VERSION,
       description = "A package providing classes to read molecules, build a hierachical tree and manipulate molecules.",
       author = 'Molecular Graphics Laboratory',
       author_email = 'sanner@scripps.edu',
       download_url = 'http://www.scripps.edu/~sanner/software/packager.html',
       url = 'http://www.scripps.edu/~sanner/software/index.html',
       packages = packages,
       py_modules = py_modules,
       data_files = data_files,
       cmdclass = {'sdist': modified_sdist,
                   'install_data': modified_install_data
                   },
       )
