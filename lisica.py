import os
import stat
import urllib2
from urllib2 import urlopen,URLError,HTTPError
import zipfile
import tkMessageBox
import json
import tempfile
import sys
import shutil
import platform
HOME_DIRECTORY=os.path.expanduser('~')
LISICA_DIRECTORY=os.path.join(HOME_DIRECTORY,"LiSiCA")




class Configuration:
    def __init__(self):
        
        self.system=platform.system()
        self.machine=platform.machine()
        self.architecture=platform.architecture()
        import struct
        self.python_bit=8 * struct.calcsize("P")
        self.python_version=platform.python_version_tuple()

    def is_os_64bit(self):
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
            print "The plugin might not be compatible with your machine"
            exe_filename="lisica"
        return exe_filename
        
    
    
class UpgraderGitlab:
 
    def __init__(self):
        self.zipFileName="archive.zip"
        self.firstVersionURL="https://git.insilab.com/insilab/lisicagui/raw/master/version.txt"        
        self.secondVersionURL="https://gitlab.com/AthiraDilip/lisicagui/raw/master/version.txt"        
        self.firstArchiveURL="https://git.insilab.com/insilab/lisicagui/repository/archive.zip"
        self.secondArchiveURL="https://gitlab.com/AthiraDilip/lisicagui/repository/archive.zip"
        self.licenseCodeGUI=""
        self.currentVersionGUI=""
        self.licenseCodeLisica=""
        self.currentVersionLisica=""
        self.latestVersionGUI=""
    
        
    def downloadInstall(self):
        try:
            print "Fetching plugin from the git server"
            urlcontent=urlopen(self.firstArchiveURL, timeout=5)
            zipcontent=urlcontent.read() 
            LiSiCAzipFile=open(self.zipFileName,'wb')
            LiSiCAzipFile.write(zipcontent)
        except:
            print "Primary git server unavailable, re-trying from secondary server"
            try:
                urlcontent=urlopen(self.secondArchiveURL, timeout=5)
                zipcontent=urlcontent.read() 
                LiSiCAzipFile=open(self.zipFileName,'wb')
                LiSiCAzipFile.write(zipcontent)   
            except HTTPError, e1:
                print "HTTP Error:", e1.code, e.read()
            except URLError, e2:
                print "URL Error:", e2.reason
        finally:
            try:
                LiSiCAzipFile.close()
            except NameError:
                pass
            try:
                urlcontent.close()
            except NameError:
                pass
        
    def extractInstall(self):
        try:
            tmp_dir= tempfile.mkdtemp();
            #with zipfile.ZipFile(zipFileName,"r") as LiSiCAzip:
            with zipfile.ZipFile(self.zipFileName,"r") as LiSiCAzip:
                for member in LiSiCAzip.namelist():
                    masterDir = os.path.dirname(member)
                    break
                LiSiCAzip.extractall(tmp_dir)
            shutil.copytree(os.path.join(tmp_dir, masterDir, "LiSiCA"),LISICA_DIRECTORY)
            shutil.rmtree(tmp_dir)
        except:
            print "installation of lisicagui failed"

    def firstUpgrade(self):
        return not os.path.isdir(LISICA_DIRECTORY)
        
    def upgrade(self):
        print "upgrading lisicagui to the latest version = ", self.latestVersionGUI
        self.downloadInstall()
        self.extractInstall()
        sys.path.append(os.path.normpath(os.path.join(LISICA_DIRECTORY,"modules")))
        import License
        License.writeToInsilabTxt(self.latestVersionGUI)
        
    def findCurrentVersion(self):
        import License
        #for GUI
        license_details=License.checkVersionGUI()
        self.licenseCodeGUI = license_details['Key']
        self.currentVersionGUI = license_details['Version']
        #for lisica program
        license_details=License.checkLicenseStatus()
        self.licenseCodeLisica = license_details['Key']
        self.currentVersionLisica = license_details['Version']

    def findLatestVersionGUI(self):
        try:
            #from insilab git server
            versionFile=urllib2.urlopen(self.firstVersionURL, timeout=5)
            self.latestVersionGUI=versionFile.read().strip()
            print "l1 = ", self.latestVersionGUI
            
        except:
            try:
                #from secondary git server
                versionFile=urllib2.urlopen(self.secondVersionURL, timeout=5)
                self.latestVersionGUI=versionFile.read().strip()
                print "l2 = ", self.latestVersionGUI
            except URLError as e:
                print e.reason
            except HTTPError as e:
                print e.read()
                print e.code
            
def run():
    
    print("Initialising LiSiCA...")
    upgraderObj = UpgraderGitlab()
    if upgraderObj.firstUpgrade():
        upgraderObj.findLatestVersionGUI()
        upgraderObj.upgrade()
        
    sys.path.append(os.path.normpath(os.path.join(LISICA_DIRECTORY,"modules")))

    configure=Configuration()
    exe_filename=configure.exe_File()
    exe_path=os.path.normpath(os.path.join(LISICA_DIRECTORY,"bin",exe_filename))
    st = os.stat(exe_path)
    os.chmod(exe_path, st.st_mode | stat.S_IEXEC)
            
    import Plugin_GUI,License

    if License.checkLicenseStatus()==None:
        active=License.activate(exe_path)
        #if not active==1:
                #pass
        
    else:  
        Plugin_GUI.main()


	

def main():
    
    run()
    
def __init__(self):
	
    self.menuBar.addmenuitem('Plugin', 'command', 'LiSiCA', label = 'LiSiCA', command=lambda s=self: main())
