#http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/252508
""" unzip.py
    Version: 1.1

    Extract a zipfile to the directoryectory provided
    It first creates the directory structure to house the files
    then it extracts the files to it.

    Sample usage:
    command line
    unzip.py -p 10 -z c:\testfile.zip -o c:\testoutput

    python class
    import unzip
    un = unzip.unzip()
    un.extract(r'c:\testfile.zip', 'c:\testoutput')
    

    By Doug Tolton
"""

import sys
import zipfile
import os
import os.path
import getopt

class unzip:
    def __init__(self, verbose = False, percent = 10):
        self.verbose = verbose
        self.percent = percent
        
    def extract(self, file, directory):
        if not directory.endswith(':') and not os.path.exists(directory):
            os.mkdir(directory)

        zf = zipfile.ZipFile(file)

        # create directoryectory structure to house files
        self._createstructure(file, directory)

        num_files = len(zf.namelist())
        percent = self.percent
        divisions = 100 / percent
        perc = int(num_files / divisions)

        # extract files to directoryectory structure
        for i, name in enumerate(zf.namelist()):
            if self.verbose == True:
                print "Extracting %s" % name
            elif perc > 0 and (i % perc) == 0 and i > 0:
                complete = int (i / perc) * percent
                print "%s%% complete" % complete

            if not name.endswith('/'):
                if os.name == 'nt': #sys.platform =='win32':
                    name1 = name.replace("/","\\")
                else:
                    name1 = name
                outfile = open(os.path.join(directory, name1), 'wb')
                outfile.write(zf.read(name))
                outfile.flush()
                outfile.close()


    def _createstructure(self, file, directory):
        self._makedirs(self._listdirs(file), directory)


    def _makedirs(self, directories, basedirectory):
        """ Create any directoryectories that don't currently exist """
        for directory in directories:
            if os.name == 'nt': #sys.platform=='win32':
                dirlist = directory.split('\\')
            else:
                dirlist = directory.split('/')
            curdirectory = basedirectory
            for direct in dirlist:
                newdir = os.path.join(curdirectory, direct)
                if not os.path.exists(newdir):
                    try:
                        os.mkdir(newdir)
                    except OSError:
                        head, tail = os.path.split(newdir)
                        os.mkdir(head)
                        os.mkdir(newdir)
                curdirectory = newdir

                
    def _listdirs(self, file):
        """ Grabs all the directoryectories in the zip structure
        This is necessary to create the structure before trying
        to extract the file to it. """
        zf = zipfile.ZipFile(file)

        directories = {}

        for name in zf.namelist():
            direct, filename = os.path.split(name)
            if os.name == 'nt': #sys.platform=='win32':
                direct = direct.replace('/', '\\')
            try:
                a = directories[direct]
            except KeyError:
                directories[direct] = 1
                
            #if sys.platform == 'win32':
            #    if name.endswith('\\'):
            #        directorys.append(name)
            #elif name.endswith('/'):
            #    directorys.append(name)
        directories = directories.keys()
        directories.sort()
        return directories

def usage():
    print """usage: unzip.py -z <zipfile> -o <targetdirectory>
    <zipfile> is the source zipfile to extract
    <targetdirectory> is the target destination

    -z zipfile to extract
    -o target location
    -p sets the percentage notification
    -v sets the extraction to verbose (overrides -p)

    long options also work:
    --verbose
    --percent=10
    --zipfile=<zipfile>
    --outdirectory=<targetdirectory>"""
    

def main():
    shortargs = 'vhp:z:o:'
    longargs = ['verbose', 'help', 'percent=', 'zipfile=', 'outdir=']

    unzipper = unzip()

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortargs, longargs)
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    zipsource = ""
    zipdest = ""

    for o, a in opts:
        if o in ("-v", "--verbose"):
            unzipper.verbose = True
        if o in ("-p", "--percent"):
            if not unzipper.verbose == True:
                unzipper.percent = int(a)
        if o in ("-z", "--zipfile"):
            zipsource = a
        if o in ("-o", "--outdir"):
            zipdest = a
        if o in ("-h", "--help"):
            usage()
            sys.exit()

    if zipsource == "" or zipdest == "":
        usage()
        sys.exit()
            
    unzipper.extract(zipsource, zipdest)

if __name__ == '__main__': main()
