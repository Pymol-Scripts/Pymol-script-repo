from Tkinter import *
from pymol import cmd
from tkFileDialog import *
def __init__(self):

    self.menuBar.addcascademenu('Plugin', 'MyPlugin', 'CASTp file selection',
                                label='CASTp pocket loader'
                                )
   
    self.menuBar.addmenuitem('MyPlugin', 'command',
                     'Remote PDB',
                       label='CASTp by PDB code',
                        command = lambda s=self : RemotePDB(s) )
                      
    self.menuBar.addmenuitem('MyPlugin', 'command', 'Get Job ID', label='CASTp by Job ID',
                             command = lambda s=self : RemoteJob(s) )
    
    self.menuBar.addmenuitem('MyPlugin', 'command',
                      'Local PDB',
                        label='CASTp from local files',
                        command = lambda s=self : LocalPDB(s) )

    self.menuBar.addcascademenu('MyPlugin', 'MyFeedback', 'FeedbackForm', label='Feedback/Bugs')
    self.menuBar.addmenuitem('MyFeedback','command', label='Leave Feedback', command = lambda s=self : Feedback(s) )
     

class Feedback:
    def __init__(self,app):
        import os
        import string
        import urllib
#        import tkCommonDiaglog
        top = Tk()
        
        F = Frame(top)
        F.pack()
        
        topLab = Label(F, text="CASTp PyMOL plug-in feedback")
        topLab.pack(padx=200)
        lHello = Label(F, text="Thank-you, your feedback is important to us!")
        lHello.pack(padx=200)
        
        A=Frame(top)
        A.pack()
        iAm = Label(A,text="I am a : ")
        aQuit = Button(A,text="Quit",command=A.quit)
        iAm.pack(side="left")
        aQuit.pack(side="left")
        
        B=Frame(top)
        B.pack()
        iUse = Label(B,text="I use the plug-in mainly for:")
        bQuit = Button(B, text="Quit", command=top.quit)
        iUse.pack(side="left")
        bQuit.pack(side="left")
        
        

#######################################################################################################
# Get pocket information from CASTp web server database.                                              #
#######################################################################################################
class RemotePDB:
    def __init__(self,app):
        import tkSimpleDialog
        import tkMessageBox
        import urllib
        import os
        import string
        remote_file = tkSimpleDialog.askstring('PDB Loader','Enter the PDB or Job ID\n\nFor PDB id\'s, you may also enter the chain.\ni.e. 1a2zA for PDB 1a2z Chain A ',parent=app.root)
        sizeof = len(remote_file)
        noerror = 1
        pdbcode = ''
        emessage = ''
        pdbfile = ''
        pocfile = ''
        infofile = ''
        jobid = ''
        
        # Gave a bad PDB code!  PDB codes must be [a-z0-9A-Z]{4}        
        if sizeof != 4 and sizeof != 5:
            tkMessageBox.showerror('Oops', remote_file + ' does not appear to be a PDB code', parent=app.root)   
            noerror = 0
        
        # Get Pocket information from CASTp web server!  Size 4 : full structure file
        # Size 5 : a single chain structure file
        if sizeof == 4:
            remote_file = remote_file.lower()
            pdbcode = remote_file
            jobid = pdbcode
            pdir = remote_file[1:3]
            path = 'http://sts.bioengr.uic.edu/castp/cast/' + pdir + '/' + remote_file + '.pdb'
            pocpath = 'http://sts.bioengr.uic.edu/castp/cast/' + pdir + '/' + remote_file + '.poc'
            infopath = 'http://sts.bioengr.uic.edu/castp/cast/' + pdir + '/' + remote_file + '.pocInfo'
        
        if sizeof == 5:
            pdbcode = remote_file[0:4]
            chident = remote_file[4:5]
            pdbcode = pdbcode.lower()
            chident = chident.upper()
            jobid = pdbcode + '.' + chident
            pdir = remote_file[1:3]
            path = 'http://sts.bioengr.uic.edu/castp/sccast/' + pdir + '/' + pdbcode + '.' + chident + '.pdb'
            pocpath = 'http://sts.bioengr.uic.edu/castp/sccast/' + pdir + '/' + pdbcode + '.' + chident + '.poc'
            infopath = 'http://sts.bioengr.uic.edu/castp/sccast/' + pdir + '/' + pdbcode + '.' + chident + '.pocInfo'
        
        # Try to retrieve the files if there are no previous errors.
        if noerror:
            pdbfile = urllib.urlretrieve(path)[0]
            pocfile = urllib.urlretrieve(pocpath)[0]
            infofile = urllib.urlretrieve(infopath)[0]
#            tkMessageBox.showerror('Gotit', pdbfile, parent=app.root)
            if(os.path.getsize(pdbfile) < 400 or os.path.getsize(pocfile) < 400 or os.path.getsize(infofile) < 400):
                emessage = pdbcode + ' is not in the CASTp database'
                noerror = 0
        
        if noerror == 0:
            tkMessageBox.showerror('Sorry', emessage, parent=app.root)            

        if noerror:
            
            # Write the contents of the pdb file to a local file. #########
            pdbin = open(pdbfile, 'r')
            pdbout = os.path.dirname(pdbfile) + os.sep + jobid + '.pdb'
            fpout = open(pdbout, 'w')
            fpout.write(pdbin.read())
            pdbin.close()
            fpout.close()
            os.remove(pdbfile)
            cmd.load(pdbout)  # Load the file
            ###############################################################

            # Read the .pocInfo file.  Create an initially empty dictionary with the
            # pocket number as the keys. Create another dictionary called pocNums
            # with the pocket number as the key and the integer pocket number as
            # the value.
            pocNums = {}
            pocDict = {}
            pocin = open(infofile, "r") 
            for line in pocin:
                stuff = line[12:16]  #Pocket Number
                stuff = stuff.strip()
                if(stuff.isdigit()):
                    pocDict[stuff] = '';        
                    idp = int(stuff)
                    pocNums[stuff] = idp

            pocin.close()
            ##############################################################
            
            # Read the .poc file.  Load the atom's of the pockets into the corresponding
            # pocket in pocDict (key = pocket number).
            pocin = open(pocfile, "r")
            for line in pocin:
                atmnum = line[6:11]
                atmnum = atmnum.strip()
                pocid = line[67:70]
                pocid = pocid.strip()
                if(atmnum.isdigit()):
                    if len(pocDict[pocid]) == 0:
                        pocDict[pocid] = atmnum
                    else:
                        pocDict[pocid] = pocDict[pocid] + "+" + atmnum

            pocin.close()
            ##############################################################################
            
            # Remove the .poc, .pocInfo and pdb file from the users temporary directory.
            os.remove(pocfile)
            os.remove(infofile)
            os.remove(pdbout)
            #############################################################################
         
            # Make an array of the pocket numbers and sort them by pocket number.
            #  Reverse the order of the sort (like CASTp webserver) such that
            #  the larger pockets will be listed first.
            pids = pocNums.values()            
            pids.sort()
            pids.reverse()
            #####################################################################
 

            
            # Load the pocket information into pyMOL!  This section is a little messy.
            #  There is a bug!  If there are too many atoms in the pocket, it can't load
            #  them all.  Need to find a way around this.
            #  Fix 1:  I reduced entries such as atom nums 1,2,3,4,5 to 1-5.  This helps,
            #          but doesn't completely fix the problem.
            counter = 0
            
            for idp in pids:  # for each pocket number (sorted)
                pid = str(idp)
                vls = pocDict[pid].split('+')
                numgrps = int(len(vls)/50)  # attempt to make groups of atoms
                stnumgrps = str(numgrps)
                currAtms = {}
                counter = 0
                
                # For each atom in the current pocket, push the atoms into an
                # array 'currAtms'.  Sort this array.
                for vl in vls:
                    ivl = int(vl)
                    currAtms[counter] = ivl
                    counter = counter + 1

                atms = currAtms.values()
                atms.sort()
                #############################################################
                
                # Here I create a new representation of the atoms.
                # If there are a group of ungapped sequential atoms
                #  I represent them as the x-y, where x is the smallest
                #  atom number in the group and y is the largest.
                #  i.e. If pocket contains 1,4,5,6,8,9,10,12..
                #       this can be represented as 1,4-6,8-10,12

                newsel = ''
                beg = ''
                currSelections = {};
                SelectionCntr = 0;
                for i in range(len(atms)):
                    if i == 0:
                        newsel = str(atms[i])
                    else:
                        newsel = newsel + "+" + str(atms[i])
        
                    if i%50 == 0 and i != 0:
                        Scntr = str(SelectionCntr)
                        tempPocket = "Pocket_" + pid + "_" + Scntr
                        currSelections[SelectionCntr] = tempPocket
                        SelectionCntr = SelectionCntr + 1
                        cmd.do("select " + tempPocket + ", id " + newsel + ",1,1")
                        newsel = ""            
            
                if newsel != "":
                    Scntr = str(SelectionCntr)
                    tempPocket = "Pocket_" + pid + "_" + Scntr
                    currSelections[SelectionCntr] = tempPocket
                    SelectionCntr = SelectionCntr + 1
                    cmd.do("select " + tempPocket + ", id " + newsel + ",0,1")
#                    cmd.do("select " + tempPocket + ", id " + newsel + ",1,1")
                    newsel = ""                    
                
                generalSelect = "select Pocket_" + pid + ", "
                
                for i in range(len(currSelections)):
                    if i==0:
                        generalSelect = generalSelect + currSelections[i]
                    else:
                        generalSelect = generalSelect + " or " + currSelections[i]
                
                generalSelect = generalSelect + ",1,1"
                cmd.do(generalSelect)
                for i in range(len(currSelections)):
                    remove = "delete " + currSelections[i]
                    cmd.do(remove)
                cmd.do("refresh")
                counter = counter + 1
                ###################################################
                
            ################################################################################################
            
############################################################################################################
# Get pocket information from the CASTp web server by job ID                                               #
############################################################################################################
class RemoteJob:
    def __init__(self,app):
        import tkSimpleDialog
        import tkMessageBox
        import urllib
        import os
        import string
        jobid = tkSimpleDialog.askstring('PDB Loader', 'Enter the Job ID given to you by the CASTp web server\nThe Job ID is case sensitive!',parent=app.root)
        pdbfile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/working/' + jobid + '.pdb')[0]
        if(os.path.getsize(pdbfile) > 400):
            pocfile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/working/' + jobid + '.poc')[0]
            pocInfofile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/working/' + jobid + '.pocInfo')[0]
           # mouthfile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/working/' + jobid + '.mouth')[0]
           # mouthInfofile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/working/' + jobid + '.mouthInfo')[0]
                                         
        else:
            os.remove(pdbfile)
            pdbfile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/uploads/' + jobid + '.pdb')[0]
            pocfile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/uploads/' + jobid + '.poc')[0]
            pocInfofile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/uploads/' + jobid + '.pocInfo')[0]
           # mouthfile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/uploads/' + jobid + '.mouth')[0]
           # mouthInfofile = urllib.urlretrieve('http://sts.bioengr.uic.edu/castp/uploads/' + jobid + '.mouthInfo')[0]

        if(os.path.getsize(pdbfile) < 400):            
            tkMessageBox.showerror('Oops!', 'Could not retrieve ' + jobid + '\nMake sure you entered it in correctly, the Job ID is case sensitive', parent=app.root)

        pdbin = open(pdbfile, 'r')
        pdbout = os.path.dirname(pdbfile) + os.sep + jobid + '.pdb'
        fpout = open(pdbout, 'w')
        fpout.write(pdbin.read())
        pdbin.close()
        fpout.close()
        os.remove(pdbfile)
        cmd.load(pdbout)
        pocNums = {}
        pocDict = {}
        pocin = open(pocInfofile, "r")
        
        for line in pocin:
            stuff = line[12:16]
            stuff = stuff.strip()
            if(stuff.isdigit()):
                pocDict[stuff] = '';        
                idp = int(stuff)
                pocNums[stuff] = idp
        pocin.close()

        pocin = open(pocfile, "r")
        for line in pocin:
            atmnum = line[6:11]
            atmnum = atmnum.strip()
            pocid = line[67:70]
            pocid = pocid.strip()
            if(atmnum.isdigit()):
                if len(pocDict[pocid]) == 0:
                    pocDict[pocid] = atmnum
                else:
                    pocDict[pocid] = pocDict[pocid] + '+' + atmnum
        pocin.close()

        os.remove(pocfile)
        os.remove(pocInfofile)
        os.remove(pdbout)
        pids = pocNums.values()
            
        pids.sort()
        pids.reverse()
        counter = 0
        for idp in pids:
            pid = str(idp)
            vls = pocDict[pid].split('+')
            currAtms = {}
            counter = 0
            for vl in vls:
                ivl = int(vl)
                currAtms[counter] = ivl
                counter = counter + 1

            atms = currAtms.values()
            atms.sort()

            newsel = ''
            beg = ''
            currSelections = {};
            SelectionCntr = 0;
            for i in range(len(atms)):
                if i == 0:
                    newsel = str(atms[i])
                else:
                    newsel = newsel + "+" + str(atms[i])

                if i%50 == 0 and i != 0:
                    Scntr = str(SelectionCntr)
                    tempPocket = "Pocket_" + pid + "_" + Scntr
                    currSelections[SelectionCntr] = tempPocket
                    SelectionCntr = SelectionCntr + 1
                    cmd.do("select " + tempPocket + ", id " + newsel + ",1,1")
                    newsel = ""            
            
            if newsel != "":
                Scntr = str(SelectionCntr)
                tempPocket = "Pocket_" + pid + "_" + Scntr
                currSelections[SelectionCntr] = tempPocket
                SelectionCntr = SelectionCntr + 1
                cmd.do("select " + tempPocket + ", id " + newsel + ",0,1")
#                cmd.do("select " + tempPocket + ", id " + newsel + ",1,1")
                newsel = ""                    
                
            generalSelect = "select Pocket_" + pid + ", "
            
            for i in range(len(currSelections)):
                if i==0:
                    generalSelect = generalSelect + currSelections[i]
                else:
                    generalSelect = generalSelect + " or " + currSelections[i]
            
            generalSelect = generalSelect + ",1,1"
            cmd.do(generalSelect)
            for i in range(len(currSelections)):
                remove = "delete " + currSelections[i]
                cmd.do(remove)
            cmd.do("refresh")
            counter = counter + 1

#######################################################################################################
# Load pocket information from files on the local machine                                             #
#######################################################################################################
class LocalPDB:
    def __init__(self,app):
        import tkMessageBox
        import tkFileDialog
        import os
        import string
#        cwd = askdirectory(title = 'Choose the directory where all of the CASTp files are located')           
        pdbfile = tkFileDialog.askopenfilename(parent=app.root, title='Open the structure file\nWithin the same directory you must have the corresponding .poc and .pocInfo files')
        wd = os.path.dirname(pdbfile) #+ os.sep + jobid + '.pdb'
        stuff = pdbfile.split('/')
        id = stuff[len(stuff)-1].split('.')
        pdbfile = wd + os.sep + id[0] + '.pdb'
        pocfile = wd + os.sep + id[0] + '.poc'
        pocInfofile = wd + os.sep + id[0] + '.pocInfo'

        cmd.load(pdbfile)
        pocNums = {}
        pocDict = {}
        pocin = open(pocInfofile, "r")
        
        for line in pocin:
            stuff = line[12:16]
            stuff = stuff.strip()
            if(stuff.isdigit()):
                pocDict[stuff] = '';        
                idp = int(stuff)
                pocNums[stuff] = idp
        pocin.close()

        pocin = open(pocfile, "r")
        for line in pocin:
            atmnum = line[6:11]
            atmnum = atmnum.strip()
            pocid = line[67:70]
            pocid = pocid.strip()
            if(atmnum.isdigit()):
                if len(pocDict[pocid]) == 0:
                    pocDict[pocid] = atmnum
                else:
                    pocDict[pocid] = pocDict[pocid] + '+' + atmnum
        pocin.close()

#        os.remove(pocfile)
#        os.remove(pocInfofile)
#        os.remove(pdbout)
        pids = pocNums.values()
           
        pids.sort()
        pids.reverse()
        counter = 0
        for idp in pids:
            pid = str(idp)
            vls = pocDict[pid].split('+')
            currAtms = {}
            counter = 0
            for vl in vls:
                ivl = int(vl)
                currAtms[counter] = ivl
                counter = counter + 1

            atms = currAtms.values()
            atms.sort()
            
            numgrps = int(len(vls)/5)  # attempt to make groups of atoms
            newsel = ''
            beg = ''
            currSelections = {};
            SelectionCntr = 0;
            for i in range(len(atms)):
                if i == 0:
                    newsel = str(atms[i])
                else:
                    newsel = newsel + "+" + str(atms[i])

                if i%50 == 0 and i != 0:
                    Scntr = str(SelectionCntr)
                    tempPocket = "Pocket_" + pid + "_" + Scntr
                    currSelections[SelectionCntr] = tempPocket
                    SelectionCntr = SelectionCntr + 1
                    cmd.do("select " + tempPocket + ", id " + newsel + ",1,1")
                    newsel = ""            
            
            if newsel != "":
                Scntr = str(SelectionCntr)
                tempPocket = "Pocket_" + pid + "_" + Scntr
                currSelections[SelectionCntr] = tempPocket
                SelectionCntr = SelectionCntr + 1
                cmd.do("select " + tempPocket + ", id " + newsel + ",0,1")
#                cmd.do("select " + tempPocket + ", id " + newsel + ",1,1")
                newsel = ""                    
                
            generalSelect = "select Pocket_" + pid + ", "
            
            for i in range(len(currSelections)):
                if i==0:
                    generalSelect = generalSelect + currSelections[i]
                else:
                    generalSelect = generalSelect + " or " + currSelections[i]
            
            generalSelect = generalSelect + ",1,1"
            cmd.do(generalSelect)
            for i in range(len(currSelections)):
                remove = "delete " + currSelections[i]
                cmd.do(remove)
            cmd.do("refresh")
            counter = counter + 1


