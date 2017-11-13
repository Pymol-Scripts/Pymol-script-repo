'''
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/lisica
 
'''

from __future__ import absolute_import
from __future__ import print_function

import os
import stat

import sys
if sys.version_info[0] < 3:
    import tkMessageBox
    import urllib2
    from urllib2 import urlopen,URLError,HTTPError
else:
    import tkinter.messagebox as tkMessageBox
    import urllib.request as urllib2
    from urllib.request import urlopen
    from urllib.error import URLError, HTTPError

import zipfile
import json
import tempfile
import sys
import shutil
import platform
from distutils.dir_util import copy_tree

HOME_DIRECTORY=os.path.expanduser('~')
LISICA_DIRECTORY=os.path.join(HOME_DIRECTORY,".lisicagui")




class Configuration:
    def __init__(self):
        
        self.system=platform.system()
        self.machine=platform.machine()
        self.architecture=platform.architecture()
        import struct
        self.python_bit=8 * struct.calcsize("P")
        self.python_version=platform.python_version_tuple()

    def is_os_64bit(self):
        if self.system == 'Windows':
            return platform.architecture()[0] == "64bit"
        else:
            return platform.machine().endswith('64')
        
    def exe_File(self):
        if self.system=='Windows':
            if self.is_os_64bit():
                exe_filename="LiSiCAx64.exe"
            else:
                exe_filename="LiSiCAx86.exe"
            
        elif self.system=='Linux':
            exe_filename="lisica"
        else:
            print("The plugin might not be compatible with your machine")
            exe_filename="lisica"
        return exe_filename
        
    
    
class UpgraderGitlab:
 
    def __init__(self):
        try:
            self.tmpDir= tempfile.mkdtemp()
            print("created temporary ", self.tmpDir)
        except:
            print("error : could not create temporary directory")
        self.zipFileName=os.path.join(self.tmpDir, "archive.zip")
        self.firstVersionURL="https://git.insilab.com/insilab/lisicagui/raw/master/.lisicagui/version.txt"        
        self.secondVersionURL="https://gitlab.com/AthiraDilip/lisicagui/raw/master/.lisicagui/version.txt"        
        self.firstArchiveURL="https://git.insilab.com/insilab/lisicagui/repository/archive.zip?ref=master"
        self.secondArchiveURL="https://gitlab.com/AthiraDilip/lisicagui/repository/archive.zip?ref=master"
        self.currentVersionGUI=""
        self.latestVersionGUI="1.0.0"


    def __del__(self):
        try:
            shutil.rmtree(self.tmpDir)
            print("deleted temporary ", self.tmpDir)
        except:
            print("error : could not remove temporary directory")
        
        
    def downloadInstall(self):
        try:
            print("Fetching plugin from the git server")
            urlcontent=urlopen(self.firstArchiveURL, timeout=5)
            zipcontent=urlcontent.read() 
            LiSiCAzipFile=open(self.zipFileName,'wb')
            LiSiCAzipFile.write(zipcontent)
        except:
            print("Primary git server unavailable, re-trying from secondary server")
            try:
                urlcontent=urlopen(self.secondArchiveURL, timeout=5)
                zipcontent=urlcontent.read() 
                LiSiCAzipFile=open(self.zipFileName,'wb')
                LiSiCAzipFile.write(zipcontent)   
            except HTTPError as e1:
                print("HTTP Error:", e1.code, e1.reason)
            except URLError as e2:
                print("URL Error:", e2.reason)
        finally:
            try:
                LiSiCAzipFile.close()
            except:
                pass
            try:
                urlcontent.close()
            except:
                pass
        
    def extractInstall(self):
        #this must be called before import LisicaGUI, otherwise
        #rmtree will fail on NFS systems due to open log file handle
        try:
            with zipfile.ZipFile(self.zipFileName,"r") as LiSiCAzip:
                for member in LiSiCAzip.namelist():
                    masterDir = os.path.dirname(member)
                    break
                LiSiCAzip.extractall(self.tmpDir)
            #copy new
            copy_tree(os.path.join(self.tmpDir, masterDir, ".lisicagui"),LISICA_DIRECTORY)
        except OSError as e:
            print("error : ", e.strerror, e.errno, e.filename)
        except shutil.Error as e:
            print("error : ", str(e))
        except:
            print("installation of lisicagui failed")

    def firstUpgrade(self):
        return not os.path.isdir(LISICA_DIRECTORY)
        
    def upgrade(self):
        print("upgrading lisicagui to the latest version = ", self.latestVersionGUI)
        self.downloadInstall()
        self.extractInstall()
        sys.path.append(os.path.normpath(os.path.join(LISICA_DIRECTORY,"modules")))
        print("Upgrade finished successfully!")
        
    def findCurrentVersion(self):
        import License
        #for GUI
        self.currentVersionGUI = License.checkVersionGUI()['Version']

    def findLatestVersionGUI(self):
        try:
            #from insilab git server
            versionFile=urllib2.urlopen(self.firstVersionURL, timeout=5)
            self.latestVersionGUI=versionFile.read().strip()
            print("l1 = ", self.latestVersionGUI)
            
        except:
            try:
                #from secondary git server
                versionFile=urllib2.urlopen(self.secondVersionURL, timeout=5)
                self.latestVersionGUI=versionFile.read().strip()
                print("l2 = ", self.latestVersionGUI)
            except HTTPError as e1:
                print("HTTP Error:", e1.code, e1.reason)
            except URLError as e2:
                print("URL Error:", e2.reason)
            
def run():
    
    print("Initialising LiSiCA...")
    try:
        sys.path.remove('')
    except:
        pass

    upgraderObj = UpgraderGitlab()

    if upgraderObj.firstUpgrade():
        upgraderObj.findLatestVersionGUI()
        upgraderObj.upgrade()
        
    sys.path.append(os.path.normpath(os.path.join(LISICA_DIRECTORY,"modules")))
    import License
    if License.checkVersionGUI()==None:
        upgraderObj.findLatestVersionGUI()
        upgraderObj.upgrade()

    del upgraderObj

    configure=Configuration()
    exe_filename=configure.exe_File()
    exe_path=os.path.normpath(os.path.join(LISICA_DIRECTORY,"bin",exe_filename))
    st = os.stat(exe_path)
    os.chmod(exe_path, st.st_mode | stat.S_IEXEC)
            

    import LisicaGUI
    LisicaGUI.main()


	

def main():
    
    run()
    
def __init__(self):
	
    self.menuBar.addmenuitem('Plugin', 'command', 'LiSiCA', label = 'LiSiCA', command=lambda s=self: main())
