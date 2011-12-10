##################################################################################
##
##Authors: Sowjanya Karnati,Michel F Sanner
##
##
################################################################################

##This command is for reporting Bug in to BugZilla Database
##
##$Id: BugReport.py,v 1.16.2.1 2011/06/06 21:06:18 sargis Exp $

import Tkinter, Pmw
from string import join
import webbrowser,os,sys
from types import StringType
from Tkinter import *
import warnings
#from htmlDoc import *


class BugReportCommand:
    
    def __init__(self,component=None):
        self.component=component
                      
    def showuploadpage_cb(self,sumcont,desccont,atcont,email_ent,
                          product="Pure Python", version="unspecified"):
        if self.component==None:
            return
       
        import urllib
        idnum = None
        ###################
        #find date and time
        ###################
        d = os.popen("date '+20%y-%m-%d'")
        date = d.readlines()[0][:-1]
        t = os.popen("date '+%H:%M:%S'")
        time = t.readlines()[0][:-1]
        deltaval = date+" "+time
        desccont += "\n sys.version: " + sys.version
        from mglutil.util.packageFilePath import findFilePath, getResourceFolderWithVersion
        registration = getResourceFolderWithVersion() + os.sep + '.registration'
        if os.path.exists(registration):
            import pickle
            try:
                desccont += "\n UserID: " + pickle.load(open(registration))['UserID'].strip()
            except Exception, inst:
                warnings.warn(inst)
                desccont += "\n Unregistered User"
        else:
            desccont += "\n Unregistered User"
        ###############################################
        ##Feeding post_bug.cgi form with required input
        ###############################################
        params = urllib.urlencode({"Bugzilla_login":"anonymous_bugzilla@yahoo.com",
                                   "Bugzilla_password":"mgltools","version":version,
                                   "rep_platform":"All","priority":"P2",
                                   "op_sys":"All","bug_severity":"normal",
                                   "bug_status":"NEW","cc":" ",
                                   "product":product,
                                   "component":"%s"  %self.component,
                                   "assigned_to":"mgltools@scripps.edu",
                                   "short_desc":"%s" %sumcont,
                                   "comment":"%s" %desccont,
                                   "bug_file_loc":" ","cc":email_ent} )
        #post
        fptr = urllib.urlopen("http://mgldev.scripps.edu/bugs/post_bug.cgi",params)
        data =fptr.readlines()
        for d in data:
            if d.endswith("</title>\n"):
                idnum = d.split(" ")[5]
                break
        if not idnum:
            return
        
        try: #this part is needed for debugging
            int(idnum)
        except:
            print data
        #for attaching files
        if len(atcont)>0:
            filelist=list(atcont)
            for f in filelist:
                if len(f)>=1:
                    params1 =  {"Bugzilla_login":"anonymous_bugzilla@yahoo.com","Bugzilla_password":"mgltools","bugid":"%i" %int(idnum),
                                "action":"insert","description":"file  attached","ispatch":"0","contenttypemethod":"autodetect","comment":" ",
                                "obsolete":"",'contenttypeselection':"",'contenttypeentry':"",'data':open(f)}
                    ###################
                    #HTTP + MULTIPART 
                    ########################
                    import urllib2
                    import MultipartPostHandler
                    opener = urllib2.build_opener(MultipartPostHandler.MultipartPostHandler)
                    urllib2.install_opener(opener)
                    req = urllib2.Request("http://mgldev.scripps.edu/bugs/attachment.cgi",params1)
                    req.add_header('Content-Type', 'text/plain')
                    response = urllib2.urlopen(req).read().strip()
         
        return idnum       
              
                
        
       
       
       
       
       
       
       
       
       ##########################################
        ###HTML PAGE
        ############################################
        #document = Tag.HTML()
        #docTitle = Tag.TITLE("Bug Submitted")
        #docHead = Tag.HEAD(docTitle)
        #docBody = [Tag.BR(),"Bug Report has been successfully submitted in to Bugzilla",Tag.BR(),Tag.BR(),"You can visit Bug at",Tag.A(href="http://mgldev.scripps.edu/bugs/show_bug.cgi?id=idnum"))
        
        #document.append([docHead,docBody])
        #fptr =open("./Pmv/BugReport.html","w")
        #HtmlElement.writeToHtml(document,HtmlFile(fptr))
        #fptr.close()
        
        #pwd = os.path.abspath("./Pmv/BugReport.html")
              
        #
        #newcont=[Tag.inputHidden("Bugzilla_login","anonymous_bugzilla@yahoo.com"),Tag.inputHidden("Bugzilla_password","mgltools"),Tag.inputHidden("version","version python2.3"),Tag.inputHidden("rep_platform","All"),Tag.inputHidden("priority","P2"),Tag.inputHidden("op_sys","All"),Tag.inputHidden("bug_severity","normal"),Tag.inputHidden("bug_status","NEW"),Tag.inputHidden("cc"," "),Tag.inputHidden("product","Pure Python"),Tag.inputHidden("component",component_n),Tag.inputHidden("assigned_to","sowjanya@scripps.edu"),Tag.inputHidden("short_desc",sumcont),Tag.inputHidden("comment",desccont),Tag.inputHidden("bug_file_loc"," ")]
        #
        ##newcont2 = [Tag.inputHidden("bugid"," "),Tag.inputHidden("data",atcont),Tag.inputHidden("description"," "),Tag.inputHidden("ispatch","1"),Tag.inputHidden("contenttypemethod","autodetect"),Tag.inputHidden("comment"," ")]
        #
        #newcont1 = [Tag.BR(),"Please press submit to submit the bug in to Bugzilla",Tag.inputSubmit(name="Submit",value="Commit")]

#newcont#=["Reporter :",Tag.TEXTAREA("Bugzilla_login",20,1,contents="anonymous_bugzilla@yahoo.com"),Tag.BR(),Tag.BR(),"Password:",Tag.inputPassword("Bugzilla_password","mgltools"),Tag.BR(),"Version: ",Tag.TEXTAREA("version",20,1,contents="version python2.3"),"          ","Platform:",Tag.TEXTAREA("rep_platform",20,1,contents="other"),Tag.BR(),"Priority:",Tag.TEXTAREA("version",20,1,contents="P2"),"          ","OS:",Tag.TEXTAREA("op_sys",20,1,contents="other"),Tag.BR(),"severity:",Tag.TEXTAREA("bug_severity",20,1,contents="normal"),Tag.BR(),"InitialState:",Tag.TEXTAREA("bug_status",20,1,contents="NEW"),Tag.BR(),"CC",Tag.TEXTAREA("cc",20,1,contents=" "),Tag.BR(),"Product:",Tag.SELECT("product",contents=[poptcont1,poptcont2],size=1),Tag.BR(),Tag.BR(),"Component:",Tag.TEXTAREA("component",20,1,contents = component_n),Tag.BR(),Tag.BR(),"AssignedTo:",Tag.TEXTAREA("assigned_to",20,1,contents = owner),Tag.BR(),"Summary:",Tag.BR(),Tag.BR(),Tag.TEXTAREA("short_desc",20,1,contents =sumcont),Tag.BR(),Tag.BR(),"Description:",Tag.BR(),Tag.TEXTAREA("comment",40,10,contents=desccont),Tag.BR(),Tag.BR(),"Attach File: ",Tag.BR(),Tag.TEXTAREA("bug_file_loc",20,1,contents = atcont),Tag.BR(),Tag.BR(),Tag.BR(),Tag.inputSubmit(name="Commit",value="Commit")]
        #

        #docBody =Tag.BODY([ Tag.heading('H2',"Report Bug"),Tag.BR(),Tag.FORM(contents=[newcont1,newcont],action="http://mgldev.scripps.edu/bugs/enter_bug.cgi",method="POST")])


        #document.append([docHead,docBody])

        ##Writing to html file
        #fptr =open("./Pmv/BugReport.html","w")
        #HtmlElement.writeToHtml(document,HtmlFile(fptr))
        #fptr.close()
        #
        #pwd = os.path.abspath("./Pmv/BugReport.html")
        ##print  sumcont,desccont
        #if len(sumcont)<=1 or len(desccont)<=1 or component_n not in ["PMV","AutoDockTools","Vision","DejaVu","FlexTree","idle","NetworkEditor","PyBabel","Pmw","ViewerFramework","Volume","tester","MolKit","mglutil","symserv"]:
        #    import tkMessageBox
        #    ok = tkMessageBox.askokcancel("Input","Please enter valid package,summary and description")
        #    return                
        #        
        #else:
        #    webbrowser.open("file://%s" %pwd)    
        
    


