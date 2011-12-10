#########################################################################
#
# Date: Jul 2003  Author: Michel Sanner, Daniel Stoffler
#
#       sanner@scripps.edu
#       stoffler@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner, Daniel Stoffler and TSRI
#
#########################################################################
#$Header: /opt/cvs/python/packages/share1.5/mglutil/util/packageFilePath.py,v 1.54.4.1 2011/06/06 17:19:38 sargis Exp $
#
#$Id: packageFilePath.py,v 1.54.4.1 2011/06/06 17:19:38 sargis Exp $

import os, string, sys, re
import warnings
import user

resourceFolder = user.home + os.sep + ".mgltools"

def getBinary(name, package):
    """Find an executable program in a Python package.  .exe extension is
assumed under windows
"""
    if (os.name == 'nt') \
      and (name [-4:]!='.exe'):
        name = name+'.exe'
    pgmfullpath = findFilePath(name, package)
    # assert it is executable
    pgmfullpath = os.path.abspath(pgmfullpath)
    if not os.access(pgmfullpath, os.X_OK):
        return None
    else:
        return pgmfullpath


def findFilePath(fileName, packageName):
    """ Find the path to the file from the package packageName"""
    mod = __import__(packageName)
    components = string.split(packageName, '.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    fullName = os.path.join(mod.__path__[0], fileName)
    if os.path.exists(fullName):
        return os.path.abspath(fullName)
    else:
        return None


def findResourceFile(module, resourceFile=None):
    """we look for the file specified in argument in the
    following directories:
    1 - current directory
    2 - user's home directory
    3 - the package to which this instance belongs to
    Returns a list of three tuples ('home',path),('currentdir',path)
    ('package',path)
    """

    if resourceFile is None:
        return
    elif resourceFile == '_visionrc':
        resourceFileLocation = {}
        resourceFileLocation['package'] = None
        resourceFileLocation['currentdir'] = None
        lResourceFolder = getResourceFolderWithVersion()
        if lResourceFolder is None:
            resourceFileLocation['home'] = None
        else:
            resourceFileLocation['home'] = getResourceFolderWithVersion() + \
                            os.sep + 'Vision' + os.sep + '_visionrc'
            if os.path.isfile(resourceFileLocation['home']) is False:
                resourceFileLocation['home'] = None            
        return resourceFileLocation
    
    resourceFileLocation = {}
    currentfile = os.path.join('.', resourceFile)
    if os.path.exists(currentfile):
        resourceFileLocation['currentdir'] = currentfile
    else:
        resourceFileLocation['currentdir'] = None

    if 'HOME' in os.environ.keys():
        home = os.environ['HOME']
        homefile = os.path.join(home, resourceFile)
        if os.path.exists(homefile):
            resourceFileLocation['home'] =  homefile
        else:
            resourceFileLocation['home'] = None
    else:
        resourceFileLocation['home'] = None

    path = __import__(module.__class__.__module__).__path__
    path = path[0]
    if path:
        packagefile = os.path.join(path, resourceFile)
        if os.path.exists(packagefile):
            resourceFileLocation['package'] =  packagefile
        else:
            resourceFileLocation['package'] = None
    else:
        resourceFileLocation['package'] = None

    return resourceFileLocation



def findAllPackages():
    """Returns a list of package names and name with path found in sys.path
"""
    packages = {}
    for p in ['.']+sys.path:
        flagline = []
        if not os.path.exists(p):
            continue
        try:
            files = os.listdir(p)
        except:
            continue
        for f in files:
            pdir = os.path.join(p, f)
            if not os.path.isdir(pdir):
                continue
            if os.path.exists( os.path.join( pdir, '__init__.py')) :
                if not packages.has_key(f):
                    packages[f] = pdir
            elif os.path.exists( os.path.join( p, '.pth') ):
                    if not packages.has_key(f):
                        packages[f] = pdir
    
    return packages            
    

def findModulesInPackage(package, name, fileNameFilters=[]):
    """
    Returns a dictionnary where the key is the path to the package or
subpackage. The value is the list of modules in which the string 'name'
was found.  Name can be a regular expression.Using '^' as a first symbol
to match string at the begining of the lines is faster.
    """

    if name[0]=='^':
        candidates = {}
        for root, dirs, files in os.walk(package):
            # remove directories not to visit
            for rem in ['CVS', 'regression', 'Tutorial', 'test', 'Doc', 'doc', 'Icons','Tests']:
                if rem in dirs:
                    dirs.remove(rem)
            # look for files that contain the string NodeLibrary
            newfiles = []
            for fi in files:
                if fi[-3:]=='.py' and not fi[0] in ['#', '.']:
                    for i in fileNameFilters:
                        if i in fi :
                            continue
                    Lines =[]        
                    f = open( os.path.join(root, fi) )
                    data = f.readlines()
                    f.close()
                    found = 0
                    Lines =filter(lambda x:x.startswith(name[1:]),data)
                    if Lines!=[]:
                        if not candidates.has_key(root): candidates[root] = []
                        candidates[root].append(fi)    
    else:  # use re
        import re
        pat = re.compile(name)
        
        candidates = {}
        for root, dirs, files in os.walk(package):
            # remove directories not to visit
            for rem in ['CVS', 'regression', 'Tutorial', 'test', 'Doc', 'doc', 'Icons','Tests']:
                if rem in dirs:
                    dirs.remove(rem)
            # look for files that contain the string NodeLibrary
            newfiles = []
            for fi in files:
                if fi[-3:]=='.py' and not fi[0] in ['#', '.']:
                    for i in fileNameFilters:
                        if i in fi :
                            continue
                    Lines =[]        
                    f = open( os.path.join(root, fi) )
                    data = f.readlines()
                    f.close()
                    found = 0
                    for line in data:
                        match = pat.search(line)
                        if match:
                            if not candidates.has_key(root): candidates[root] = []
                            candidates[root].append(fi)
                            break
    return candidates
    

def findModulesInDirectory(directory, name):
    """Returns a list of modules for a given directory, in which the string
'name' was found."""
     
    #pat = re.compile(name)
    candidates = {}
    olddir = os.getcwd()
    os.chdir(directory)
    from glob import glob
    files = glob("*.py")
    # look for files that contain the string NodeLibrary
    for fi in files:
        Lines =[]
        f = open( fi )
        data = f.readlines()
        f.close()
        found = 0
        Lines =filter(lambda x:x.startswith(name),data)
        if Lines!=[]:
            if not candidates.has_key(root): candidates[root] = []
            candidates[root].append(fi)    
            
    os.chdir(olddir)
    return candidates
    

def getObjectFromFile(absfile, obj):
    """import the object obj from absfile and return it"""
    direct, file = os.path.split(absfile)
    sys.path.insert(0, direct)
    modulepath = file[:-3].replace(os.path.sep, '.')
    mod = __import__(modulepath, globals(), locals(), [obj])
    components = string.split(modulepath, '.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    sys.path = sys.path[1:]
    return getattr(mod, obj)


def setResourceFolder(folderName):
    """ Set the name of the resourceFolder
"""
    global resourceFolder
    resourceFolder = folderName
    if os.path.isdir(resourceFolder) is False:
        try:
            os.mkdir(resourceFolder)
        except Exception, e:
            txt = "Resource Folder can't be created in home directory, now using tmp directory %s %s" %(resourceFolder, e)
            warnings.warn(txt)
            from tempfile import gettempdir
            resourceFolder = gettempdir() + os.path.sep + ".mgltools"
            if not os.path.isdir(resourceFolder):
                try:
                    os.mkdir(resourceFolder)
                except:
                    txt = "Running without Resource Folder as it can't be created %s" %resourceFolder
                    warnings.warn(txt)
                    resourceFolder = None
    return resourceFolder


def getResourceFolder():
    """ Returns MGLTools resource folder,
create it if necessary, 
returns None if it doesn't exist and can't be created
"""
    global resourceFolder
    if resourceFolder is None:
        return None
    if os.path.isdir(resourceFolder) is False:
        setResourceFolder(resourceFolder)
    return resourceFolder


def getResourceFolderWithVersion():
    """ Returns MGLTools resource folder,
create it if necessary, 
returns None if it doesn't exist and can't be created
"""
    old_rc = getResourceFolder()
    if old_rc is None:
        return None
    
    from Support.version import __version__
    Version = __version__

    lRessourceFolder = old_rc + os.sep + Version

    if os.path.isdir(lRessourceFolder):
        return lRessourceFolder
    else:
        try:
            os.mkdir(lRessourceFolder)                    
            return lRessourceFolder
        except Exception, inst:
            print inst
            print "Cannot create the Resource Folder %s" %lRessourceFolder
            import tempfile
            tmpDir = tempfile.gettempdir()
            print "Using %s insterad" %tmpDir
            return tmpDir
        
def which (filename):
    """Source http://mail.python.org/pipermail/python-list/2002-August/157829.html"""
    if os.access(filename, os.X_OK):
            return filename    
    if not os.environ.has_key('PATH') or os.environ['PATH'] == '':
        p = os.defpath
    else:
        p = os.environ['PATH']

    pathlist = p.split (os.pathsep)

    for path in pathlist:
        f = os.path.join(path, filename)
        if os.access(f, os.X_OK):
            return f
        fAlternative = f + '.exe'
        if os.access(fAlternative, os.X_OK):
            return fAlternative
        fAlternative = f + '.sh'
        if os.access(fAlternative, os.X_OK):
            return fAlternative
        fAlternative = f + '.bat'
        if os.access(fAlternative, os.X_OK):
            return fAlternative        
    return None