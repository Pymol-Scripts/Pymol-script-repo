import Tkinter
from Tkinter import *
import os
import subprocess
import tkMessageBox



HOME_DIRECTORY=os.path.expanduser('~')
LISICA_DIRECTORY=os.path.join(HOME_DIRECTORY,".lisicagui")
versionFile = os.path.join(LISICA_DIRECTORY,"version.txt")

def writeToInsilabTxt(latestVersion):
    # needed for backward compatibility
    pass

def checkLicenseStatus():
    # needed for backward compatibility
    return {'Key':"FREE_FOR_ACADEMIC_USE",'Version':"1.0.0"}

def checkVersionGUI():
    if os.path.isfile(versionFile):
        with open(versionFile) as vFile:
            lines = vFile.readlines()
        return {'Key':"FREE_FOR_ACADEMIC_USE",'Version':lines[0].strip()}
    return None

    
