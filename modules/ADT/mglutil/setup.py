#!/usr/bin/env python

from distutils.core import setup
from distutils.command.sdist import sdist
from distutils.command.install_data import install_data
from glob import glob
import os

########################################################################
# Had to overwrite the prunefile_list method of sdist to not
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
packages = ['mglutil', 'mglutil.TestUtil', 'mglutil.Tests', 'mglutil.gui',
            'mglutil.gui.BasicWidgets', 'mglutil.gui.BasicWidgets.Tk',
            'mglutil.gui.BasicWidgets.Tk.icons',
            'mglutil.gui.BasicWidgets.Tk.Tests',
            'mglutil.gui.BasicWidgets.Tk.TreeWidget',
            'mglutil.gui.BasicWidgets.Tk.TreeWidget.icons',
            'mglutil.gui.BasicWidgets.Tk.TreeWidget.Tests',
            'mglutil.gui.InputForm', 'mglutil.gui.InputForm.Tk',
            'mglutil.gui.InputForm.Tk.Tests',
            'mglutil.gui.Misc','mglutil.gui.Misc.Tk',
            'mglutil.math', 'mglutil.util',
            'mglutil.regression',
            'mglutil.TestUtil.Tests',
            'mglutil.TestUtil.Tests.Data',
            'mglutil.web',
            'mglutil.web.regression',
            'mglutil.web.services',
            'mglutil.util.Tests',
            'mglutil.web.Tests',
            'mglutil.math.Tests',
            'mglutil/TestUtil/bin',
            'mglutil/splashregister',
            'mglutil/gui/BasicWidgets/Tk/trees',
            'mglutil/gui/Spline', 'mglutil/gui/Spline/Tests']

# list of the python modules which are not part of a package
py_modules = []

# list of the files that are not python packages but are included in the
# distribution and need to be installed at the proper place  by distutils.
# The list in MANIFEST.in lists is needed for including those files in
# the distribution, data_files in setup.py is needed to install them
# at the right place.


data_files = []

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

os.path.walk("mglutil", getDataFiles, data_files)
#print "data_files:", data_files, len(data_files)

from version import VERSION
setup (name = "mglutil",
       version = VERSION,
       description = "Molecular Graphics Laboratory utility collection",
       author = "Molecular Graphics Laboratory",
       author_email = "sanner@scripps.edu",
       download_url = 'http://www.scripps.edu/~sanner/software/packager.html',
       url = "http://www.scripps.edu/~sanner/software/index.html",
       license = "to be specified",
       packages = packages,
       cmdclass = {'install_data': modified_install_data,
                   'sdist': modified_sdist,},
       data_files = data_files
       )
