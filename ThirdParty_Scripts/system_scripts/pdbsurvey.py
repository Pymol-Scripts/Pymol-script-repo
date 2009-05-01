from pymol import cmd
 
def pdbsurvey(days=50):
 
    """USAGE : pdbsurvey (<days>)
    Surveys the updates added to the PDB (ftp.rcsb.org) in the last
    50 days (or otherwise specified when calling this function) for
    entries that contain the words specified in the file
    keywords.txt.
    """
    print days
 
    import ftplib
    import time
    import os
    import string
 
 
 
    def todaymerge():
        """Puts today's date in a pdb format string.
        """
        date=time.localtime()
        fyear="%i" %(date[0])
        fmonth="%i" %(date[1])
        if date[1]<10:
            fmonth="0"+"%i" %(date[1])
        fday="%i" %(date[2])
        if date[2]<10:
            fday="0"+"%i" %(date[2])
        dateS=fyear+fmonth+fday
        return dateS
 
    def file2list(filename):
        """Low-level routine to brainlessly implement
        file.read().
        """
        fq=open(filename,'rb')
        linesS=fq.read()
        fq.close()
        LIST=linesS.splitlines()
        return LIST
 
    def connect2pdb():
        """Opens an anonymous socket to ftp://ftp.rcsb.org
        """
        f=ftplib.FTP()
        f.connect ('ftp.rcsb.org')
        f.login ()
        print "Remote connection established","\n"
        return f
 
    def decrementdate(dateS):
        """given a string date (pdb format yyyymmdd)
        this routine returns a string of the day before
        (sadly assuming that every month has 31 days, but
        no big deal here...).
        """
        #decompose dateS into components
        yearS=dateS[0]+dateS[1]+dateS[2]+dateS[3]
        monthS=dateS[4]+dateS[5]
        dayS=dateS[6]+dateS[7]
 
        #convert each into integers
        yearI=int(yearS)
        monthI=int(monthS)
        dayI=int(dayS)
 
        #actual stuff
        dayI=dayI-1
        if dayI==0:
            dayI=31
            monthI=monthI-1
            if monthI==0:
                monthI=12
                yearI=yearI-1
        dayS="%i" %(dayI)
        monthS="%i" %(monthI)
        yearS="%i" %(yearI)
        if dayI<10:
            dayS="0"+dayS
        if monthI<10:
            monthS="0"+monthS
        #and finally...
        dateS=yearS+monthS+dayS
        return dateS
 
    def findlastdir(dateS,f,days):
        """Puts the names of the "recent" directories in the
        list named "directoriesL".
        """
        directoriesL=['']
        j=p=0
        while p<days:
            dateS=decrementdate(dateS)
            attempt="/pub/pdb/data/status/"+dateS
            try :
                f.cwd(attempt)
                directoriesL[j:j]=[attempt]
                j=j+1
            except :
                pass
            p=p+1
        directoriesL.pop()
        return directoriesL
 
    def compilinfile(directoriesL,f):
        """lists all structures in the added.pdb files
        contained in the directories specified in directoriesL
        """
        command="RETR added.pdb"
        handle=open("donotedit.dat","wrb")
        for k in directoriesL:
            f.cwd(k)
            print "Currently in directory ",f.pwd()
            f.retrbinary(command,handle.write)
        handle.close()
        return len(directoriesL)
 
    def listparser():
        """Extracts the pdbids from donotedit.dat file,
        and stacks them into the list pdbidsL
        """
        linesL=file2list("donotedit.dat")
        pdbidsL=[]
        for iter in linesL:
            pdbidsL.append(iter[57:61])
        for iter in pdbidsL:
            iter=string.lower(iter)
        pdbidsL.sort()
        return pdbidsL
 
    def currentrelease(f):
        """Stores the content of cmpd_res.idx file
        This file contains the equivalencies pdbid<->title
        for all current entries of the PDB.
        """
        command="RETR cmpd_res.idx"
        f.cwd("/pub/pdb/derived_data/index/")
        print "Currently in directory ",f.pwd()
        fq=open("dictionnary.dat",'wrb')
        f.retrbinary(command,fq.write)
        fq.close()
        dictL=file2list("dictionnary.dat")
        return dictL
 
    def extract(pdbidsL,dictL):
        """Populates dictionnaryD with pdb entries found in the
        latest releases.
        """
        dictionnaryD={}
        problemL=[]
        extractL=[dictionnaryD,problemL]
        for i in dictL:
            tempS=string.lower(i[0:4])
            for ii in pdbidsL:
                if ii == tempS:
                    title=i[14:216]
                    extractL[0][ii]=title
        if len(extractL[0].keys()) != len(pdbidsL):
            print "Dimension mismatch, seeking troublemaker..."
            for i in pdbidsL:
                equiv=0
                for ii in extractL[0].keys():
                    if i==ii:
                        equiv=equiv+1
                if equiv==0:
                    extractL[1].append(i)
        return extractL
 
    def disconnectpdb(f):
        """Diconnects the current ftp session
        """
        f.quit()
        print "Remote connection terminated","\n"
        return f
 
    def releventries(dictionnaryD):
        """Generates a cleaned dictionnary with only entries
        that have one or more keywords specified in the local
        user-defined keywords.txt file
        """
        keywL=file2list("keywords.txt")
        relevdicD={}
        for i in keywL:
            for elem in dictionnaryD.keys():
                temp=dictionnaryD[elem]
                if temp.find(i) != -1:
                    relevdicD[elem]=temp
        return relevdicD
 
    def diskcleanup(filelist=["donotedit.dat","dictionnary.dat"]):
        """Lo-level disk cleanup to free up memory without the user
        """
        for filename in filelist:
            command='DEL '+filename
            os.system(command)
        return "clean"
 
 
 
 
    print "Welcome in the auto-PDB updater !"
 
    print "Survey of updates made since",days,"days ago."
 
    print "Acquisition of local time..."
    dateS=todaymerge()                                                 #Initializes dateS
    print "today is ",dateS
    print "Connecting to remote ftp server..."
    f=connect2pdb()                                                    #Connect anonymously to ftp.rcsb.org
 
    print "Acquisition of latest added remote directories..."
    directoriesL=findlastdir(dateS,f,days)                             #Lists recent directories in directoriesL
    if len(directoriesL)==0:
        print "No updates have been found since",days,"ago. Starting over with 50 days ago."
        directoriesL=findlastdir(dateS,f,50)
 
    print "Acquisition of latest addedremote files..."
    updatesnumberI=compilinfile(directoriesL,f)                        #Concatenates the corresponding added.pdb into donotedit.dat
 
    print "Parsing of latest entries..."
    pdbidsL=listparser()                                               #Recent names now present in the pdbidsL list (one name per element)
 
    print "Acquisition of the current pdb distribution..."
    dictL=currentrelease(f)                                            #Populates dictL with the current entries of the PDB
 
    print "Parsing of the current pdb distribution into [code,title] tuples..."
    extractL=extract(pdbidsL,dictL)                                    #generates the dictionnary of latest releases key:PDBid ; definition:pdbtitle
 
    print "Disconnection from the remote ftp server..."
    f=disconnectpdb(f)                                                 #Closes the ftp instance
 
    print "Extraction of the relevant entries..."
    relevdicD=releventries(extractL[0])                               #Generates a subset of dictionnary D with criterion being "has keywords contained in keywords.txt in its title"
 
    print "Cleaning program-generated temporary files..."
    clean=diskcleanup()                                                #Cleans the mess generated by the program
 
    reportL=[]
    reportL.append("\n")
    reportL.append("###############REPORT########################################\n")
    reportL.append("\n")
    lendictS="%i" %(len(dictL))
    chmilblik = 'The current pdb version (as of '+dateS+") has "+lendictS+" entries.\n"
    reportL.append(chmilblik)
    line="The most recent directory is : "+directoriesL[0]+".\n"
    reportL.append(line)
    updatesnumberS="%i" %(updatesnumberI)
    entriesnumber="%i" %(len(extractL[0].keys()))
    line="The "+updatesnumberS+" last updates ("+entriesnumber+" entries) have been examined.\n"
    reportL.append(line)
    diclengthS="%i" %(len(relevdicD.keys()))
    line=diclengthS+" are relevant to you :\n"
    reportL.append(line)
    for i in relevdicD.keys():
        entry=i+" : "+relevdicD[i]+"\n"
        reportL.append(entry)
    problemS=""
    for i in extractL[1]:
        problemS=i+";"+problemS
    problemS="["+problemS
    problemS=problemS.strip(";")
    problemS=problemS+"]"
    lineS="The entries "+problemS+" raised problems,"
    reportL.append(lineS)
    reportL.append("they should be examined manually.")
    reportL.append("\n")
    reportL.append("###############END OF REPORT#################################\n")
    report=open("report.aut","w")
    for elem in reportL:
        print elem
        elem=elem+'\n'
        report.writelines(elem)
    report.close()
    command2='start keywords.txt'
    command3='start report.aut'
    os.system(command2)
    os.system(command3)
 
cmd.extend("pdbsurvey",pdbsurvey)
