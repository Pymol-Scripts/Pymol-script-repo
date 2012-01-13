"""
    Driver for PDB2PQR

    This module takes a PDB file as input and performs optimizations
    before yielding a new PDB-style file as output.

    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

    Parsing utilities provided by Nathan A. Baker (Nathan.Baker@pnl.gov)
    Pacific Northwest National Laboratory

    Copyright (c) 2002-2011, Jens Erik Nielsen, University College Dublin; 
    Nathan A. Baker, Battelle Memorial Institute, Developed at the Pacific 
    Northwest National Laboratory, operated by Battelle Memorial Institute, 
    Pacific Northwest Division for the U.S. Department Energy.; 
    Paul Czodrowski & Gerhard Klebe, University of Marburg.

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
        * Neither the names of University College Dublin, Battelle Memorial Institute,
          Pacific Northwest National Laboratory, US Department of Energy, or University
          of Marburg nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

"""

__date__  = "5 April 2010"
__author__ = "Todd Dolinsky, Nathan Baker, Jens Nielsen, Paul Czodrowski, Jan Jensen, Samir Unni, Yong Huang"
__version__ = "1.7.1"


import string
import sys
import getopt
import os
import time
import httplib
import re
import glob
import tempfile
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
#from src import server
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *
from src.routines import *
from src.protein import *
from src.server import *
from src.hydrogens import *
from src.aconf import *
from StringIO import *
from main import *

def printHeader(pagetitle,have_opal=None,jobid=None):
    """
        Function to print html headers
    """
    if jobid:
        if have_opal:
            print "Location: querystatus.cgi?jobid=%s&typeofjob=opal\n" % (jobid,typeOfJob)
        else:
            print "Location: querystatus.cgi?jobid=%s&typeofjob=local\n" & (jobid,typeOfJob)

    #print "Content-type: text/html\n"
    print "<HTML>"
    print "<HEAD>"
    print "\t<TITLE>%s</TITLE>" % pagetitle
    print "\t<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET
    print "</HEAD>"
    return

def redirector(name):
    """
        Prints a page which redirects the user to querystatus.cgi and writes starting time to file
    """

    utilities.appendToLogFile(name, 'pdb2pqr_start_time', str(time.time()))

    string = ""
    string+= "<html>\n"
    string+= "\t<head>\n"
    string+= "\t\t<meta http-equiv=\"Refresh\" content=\"0; url=%squerystatus.cgi?jobid=%s&calctype=pdb2pqr\">\n" % (WEBSITE, name)
    string+= "\t</head>\n"
    string+= "</html>\n"
    return string

def sanitizeFileName(fileName):
    fileName = re.split(r'[/\\]',fileName)[-1]
    fileName = fileName.replace(' ', '_')
    fileName = fileName.replace('-', '_')
    return fileName

class WebOptionsError(Exception):
    pass

class WebOptions(object):
    '''Helper class for gathering and querying options selected by the user'''
    def __init__(self, form):
        '''Gleans all information about the user selected options and uploaded files.
        Also validates the user input. raises WebOptionsError if there is any problems.'''
        
        #options to pass to runPDB2PQR
        self.runoptions = {}
        #Additional options to pass to google analytics along with the run options.
        #These are included in has_key(), __contains__(), and __getitem__() calls.
        self.otheroptions = {}
        
        self.runoptions['debump'] = form.has_key("DEBUMP")
        self.runoptions['opt'] = form.has_key("OPT")
        
        if form.has_key('FF'):
            self.ff = form["FF"].value.lower()
        else:
            raise WebOptionsError('Force field type missing from form.')
        
        if form.has_key("PDBID") and form["PDBID"].value and form["PDBSOURCE"].value == 'ID':
            self.pdbfile = getPDBFile(form["PDBID"].value)
            if self.pdbfile is None:
                raise WebOptionsError('The pdb ID provided is invalid.')
            self.pdbfilestring = self.pdbfile.read()
            self.pdbfile = StringIO(self.pdbfilestring)
            self.pdbfilename = form["PDBID"].value
        elif form.has_key("PDB") and form["PDB"].filename and form["PDBSOURCE"].value == 'UPLOAD':
            self.pdbfilestring = form["PDB"].value
            self.pdbfile = StringIO(self.pdbfilestring)
            self.pdbfilename = sanitizeFileName(form["PDB"].filename)
        else:
            raise WebOptionsError('You need to specify a pdb ID or upload a pdb file.')
            
        if form.has_key("PROPKA"):
            if not form.has_key('PH'):
                raise WebOptionsError('Please provide a pH value.')
            
            phHelp = 'Please choose a pH between 0.0 and 14.0.'
            try:
                ph = float(form["PH"].value)
            except ValueError:
                raise WebOptionsError('The pH value provided must be a number!  ' + phHelp)
            if ph < 0.0 or ph > 14.0: 
                text = "The entered pH of %.2f is invalid!  " % ph
                text += phHelp
                raise WebOptionsError(text)
            self.runoptions['ph'] = ph
                 
        self.otheroptions['apbs'] = form.has_key("INPUT")
        self.otheroptions['whitespace'] = form.has_key("WHITESPACE")
        
        if self.ff == 'user':
            if form.has_key("USERFF") and form["USERFF"].filename:
                self.userfffilename = sanitizeFileName(form["USERFF"].filename)
                self.userffstring = form["USERFF"].value
                self.runoptions['userff'] = StringIO(form["USERFF"].value)
            else:
                text = "A force field file must be provided if using a user created force field."
                raise WebOptionsError(text)
                
            if form.has_key("USERNAMES") and form["USERNAMES"].filename:
                self.usernamesfilename = sanitizeFileName(form["USERNAMES"].filename)
                self.usernamesstring = form["USERNAMES"].value
                self.runoptions['usernames'] = StringIO(form["USERNAMES"].value)
            else:
                text = "A names file must be provided if using a user created force field."
                raise WebOptionsError(text)
            
        if form.has_key("FFOUT") and form["FFOUT"].value != "internal":
            self.runoptions['ffout'] = form["FFOUT"].value
                
        self.runoptions['chain'] = form.has_key("CHAIN")
        self.runoptions['typemap'] = form.has_key("TYPEMAP")
        self.runoptions['neutraln'] = form.has_key("NEUTRALN")
        self.runoptions['neutralc'] = form.has_key("NEUTRALC")
        
        if (self.runoptions['neutraln'] or self.runoptions['neutraln']) and \
            self.ff != 'parse':
            raise WebOptionsError('Neutral N-terminus and C-terminus require the PARSE forcefield.')
        
        if form.has_key("LIGAND") and form['LIGAND'].filename:
            self.ligandfilename=sanitizeFileName(form["LIGAND"].filename)
            ligandfilestring = form["LIGAND"].value
            # for Windows-style newline compatibility for pdb2pka
            ligandfilestring = ligandfilestring.replace('\r\n', '\n')
            self.ligandfilestring = ligandfilestring.replace('\r', '\n')
            
            self.runoptions['ligand'] = StringIO(self.ligandfilestring)
            
        if self.pdbfilename[-4:]==".pdb":
            self.pqrfilename = "%s.pqr" % self.pdbfilename[:-4]
        else:
            self.pqrfilename = "%s.pqr" % self.pdbfilename
            
        self.runoptions['verbose'] = False
        
    def getLoggingList(self):
        '''Returns a list of options the user has turned on.
        Used for logging jobs later in usage.txt'''
        results = []
        
        for key in self:
            if self[key]:
                results.append(key)
                
        return results
        
    def getRunArguments(self):
        '''Returns argument suitable for runPDB2PQR'''
        return self.runoptions.copy()
    
    def getOptions(self):
        '''Returns all options for reporting to google analytics and local job logging.
        For reporting purposes only!'''
        options = self.runoptions.copy()
        options.update(self.otheroptions)
        
        #we need stuff that is easily converted to and from a string.
        #This is for reporting only.
        if 'ligand' in self.runoptions:
            options['ligand'] = self.ligandfilename
            
        if 'usernames' in self.runoptions:
            options['usernames'] = self.usernamesfilename
            
        if 'userff' in self.runoptions:
            options['userff'] = self.userfffilename
        
        return options
    
    def getCommandLine(self):
        commandLine = []
        
        if not self.runoptions['debump']:
            commandLine.append('--nodebump')
            
        if not self.runoptions['opt']:
            commandLine.append('--noopt')
            
        if 'ph' in self.runoptions:
            commandLine.append('--with-ph=%s' % self.runoptions['ph'])
        
        if self.otheroptions['apbs']:
            commandLine.append('--apbs-input')
            
        if self.otheroptions['whitespace']:
            commandLine.append('--whitespace')
            
        if 'userff' in self.runoptions and self.ff == 'user':
            commandLine.append('--userff=%s' % self.userfffilename)
            commandLine.append('--usernames=%s' % self.usernamesfilename)
        else:
            commandLine.append('--ff=%s' % self.ff)
            
        if 'ffout' in self.runoptions:
            commandLine.append('--ffout=%s' % self.runoptions['ffout'])
            
        for o in ('chain', 'typemap', 'neutraln', 'neutralc'):
            if self.runoptions[o]:
                commandLine.append('--' + o)
            
        if 'ligand' in self.runoptions:
            commandLine.append('--ligand=%s' % self.ligandfilename)
            
        commandLine.append(self.pdbfilename)
        
        commandLine.append(self.pqrfilename)
        
        return ' '.join(commandLine)
    
    def __contains__(self, item):
        '''Helper for checking for the presence of an option'''
        return item in self.runoptions or item in self.otheroptions
    
    def has_key(self, item):
        '''Helper for checking for the presence of an option'''
        return item in self.runoptions or item in self.otheroptions
    
    def __iter__(self):
        for key in self.runoptions:
            yield key
        
        for key in self.otheroptions:
            yield key
            
    def __getitem__(self, key):
        return self.runoptions[key] if key in self.runoptions else self.otheroptions[key]

def handleOpal(weboptions):
    '''
        Handle opal based run.
    '''
    # Opal-specific import statments
    from AppService_client import AppServiceLocator, getAppMetadataRequest, launchJobRequest, launchJobBlockingRequest, getOutputAsBase64ByNameRequest
    from AppService_types import ns0
    from ZSI.TC import String

    inputFiles = []
    
    if 'userff' in weboptions:
        ffFile = ns0.InputFileType_Def('inputFile')
        ffFile._name = weboptions.userfffilename
        ffFile._contents = weboptions.userffstring
        inputFiles.append(ffFile)
        
    if 'usernames' in weboptions:
        namesFile = ns0.InputFileType_Def('inputFile')
        namesFile._name = weboptions.usernamesfilename
        namesFile._contents = weboptions.usernamesstring
        inputFiles.append(namesFile)
        
    if 'ligand' in weboptions:
        ligandFile = ns0.InputFileType_Def('inputFile')
        ligandFile._name = weboptions.ligandfilename
        ligandFile._contents = weboptions.ligandfilestring
        inputFiles.append(ligandFile)
        
    pdbOpalFile = ns0.InputFileType_Def('inputFile')
    pdbOpalFile._name = weboptions.pdbfilename
    pdbOpalFile._contents = weboptions.pdbfilestring
    inputFiles.append(pdbOpalFile)
     
    # launch job   
    appLocator = AppServiceLocator()
    appServicePort = appLocator.getAppServicePort(PDB2PQR_OPAL_URL)
    
    req = launchJobRequest()
    req._argList = weboptions.getCommandLine()
    req._inputFile=inputFiles

    try:
        resp=appServicePort.launchJob(req)
    except Exception, e:
        printHeader("PDB2PQR Job Submission - Error")
        print "<BODY>\n<P>"
        print "There was an error with your job submission<br>"
        print "</P>"
        print "<script type=\"text/javascript\">"
        print "var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");"
        print "document.write(unescape(\"%3Cscript src=\'\" + gaJsHost + \"google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E\"));"
        print "</script>"
        print "<script type=\"text/javascript\">"
        print "try {"
        print "var pageTracker = _gat._getTracker(\"UA-11026338-3\");"
        for key in weboptions:
            print "pageTracker._trackPageview(\"/main_cgi/has_%s_%s.html\");" % (key, weboptions[key])
        print "pageTracker._trackPageview();"
        print "} catch(err) {}</script>"
        print "</BODY>"
        print "</HTML>"
        sys.exit(2)
    
    try:
        starttime = time.time()
        name = setID(starttime)
        
        #Some job parameters logging.
        os.makedirs('%s%s%s' % (INSTALLDIR, TMPDIR, name))
        apbsInputFile = open('%s%s%s/apbs_input' % (INSTALLDIR, TMPDIR, name),'w')
        apbsInputFile.write(str(weboptions["apbs"]))
        apbsInputFile.close()
        
        typemapInputFile = open('%s%s%s/typemap' % (INSTALLDIR, TMPDIR, name),'w')
        typemapInputFile.write(str(weboptions["typemap"]))
        typemapInputFile.close()
        
        pdb2pqrOpalJobIDFile = open('%s%s%s/pdb2pqr_opal_job_id' % (INSTALLDIR, TMPDIR, name), 'w')
        pdb2pqrOpalJobIDFile.write(resp._jobID)
        pdb2pqrOpalJobIDFile.close()
        
        print redirector(name)
        
        # Recording CGI run information for PDB2PQR Opal
        pdb2pqrOpalLogFile = open('%s%s%s/pdb2pqr_opal_log' % (INSTALLDIR, TMPDIR, name), 'w')
        pdb2pqrOpalLogFile.write(str(weboptions.getOptions())+'\n'+
                                 str(weboptions.ff)+'\n'+
                                 str(os.environ["REMOTE_ADDR"]))
        pdb2pqrOpalLogFile.close()

    except StandardError, details:
        print details
        createError(name, details)

def handleNonOpal(weboptions):
    """
        Handle non opal run.
    """
 
    pdblist, errlist = readPDB(weboptions.pdbfile)
    
    dummydef = Definition()
    dummyprot = Protein(pdblist, dummydef)
    if len(pdblist) == 0 and len(errlist) == 0:
        text = "Unable to find PDB file - Please make sure this is "
        text += "a valid PDB file ID!"
        #print "Content-type: text/html\n"
        print text
        sys.exit(2)
    elif dummyprot.numAtoms() > MAXATOMS and weboptions["opt"] == True:
        text = "<HTML><HEAD>"
        text += "<TITLE>PDB2PQR Error</title>"
        text += "<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">" % STYLESHEET
        text += "</HEAD><BODY><H2>PDB2PQR Error</H2><P>"
        text += "Due to server limits, we are currently unable to optimize "
        text += "proteins of greater than MAXATOMS atoms on the server (PDB2PQR "
        text += "found %s atoms in the selected PDB file).  If you " % dummyprot.numAtoms()
        text += "want to forgo optimization please try the server again.<P>"
        text += "Otherwise you may use the standalone version of PDB2PQR that "
        text += "is available from the <a href=\"http://pdb2pqr.sourceforge.net\">"
        text += "PDB2PQR SourceForge project page</a>."
        text += "<script type=\"text/javascript\">"
        text += "var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");"
        text += "document.write(unescape(\"%3Cscript src=\'\" + gaJsHost + \"google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E\"));"
        text += "</script>"
        text += "<script type=\"text/javascript\">"
        text += "try {"
        text += "var pageTracker = _gat._getTracker(\"UA-11026338-3\");"
        for key in weboptions:
            text += "pageTracker._trackPageview(\"/main_cgi/has_%s_%s.html\");" % (key, weboptions[key])
        text += "pageTracker._trackPageview();"
        text += "} catch(err) {}</script>"
        text += "</BODY></HTML>"
        #print "Content-type: text/html\n"
        print text
        sys.exit(2)

    try:
        starttime = time.time()
        name = setID(starttime)

        #Some job parameters logging.
        os.makedirs('%s%s%s' % (INSTALLDIR, TMPDIR, name))
        apbsInputFile = open('%s%s%s/apbs_input' % (INSTALLDIR, TMPDIR, name),'w')
        apbsInputFile.write(str(weboptions["apbs"]))
        apbsInputFile.close()
        
        typemapInputFile = open('%s%s%s/typemap' % (INSTALLDIR, TMPDIR, name),'w')
        typemapInputFile.write(str(weboptions["typemap"]))
        typemapInputFile.close()

        statusfile = open('%s%s%s/pdb2pqr_status' % (INSTALLDIR, TMPDIR, name), 'w')
        statusfile.write('running')
        statusfile.close()


        pid = os.fork()
        if pid:
            print redirector(name)
            sys.exit()
        else:
            currentdir = os.getcwd()
            os.chdir("/")
            os.setsid()
            os.umask(0)
            os.chdir(currentdir)
            os.close(1) # not sure if these
            os.close(2) # two lines are necessary
            
            pqrpath = '%s%s%s/%s.pqr' % (INSTALLDIR, TMPDIR, name, name)
            
            orig_stdout = sys.stdout
            orig_stderr = sys.stderr
            sys.stdout = open('%s%s%s/pdb2pqr_stdout.txt' % (INSTALLDIR, TMPDIR, name), 'w')
            sys.stderr = open('%s%s%s/pdb2pqr_stderr.txt' % (INSTALLDIR, TMPDIR, name), 'w')
            
            header, lines, missedligands = runPDB2PQR(pdblist, 
                                                      weboptions.ff,
                                                      outname = pqrpath,
                                                      **weboptions.getRunArguments())
            
            sys.stdout.close()
            sys.stderr.close()
            sys.stdout = orig_stdout
            sys.stderr = orig_stderr

            endtimefile = open('%s%s%s/pdb2pqr_end_time' % (INSTALLDIR, TMPDIR, name), 'w')
            endtimefile.write(str(time.time()))
            endtimefile.close()

            pqrfile = open(pqrpath, "w")
            pqrfile.write(header)
            
            whitespace = weboptions.otheroptions['whitespace']
            for line in lines:
                # Adding whitespaces if --whitespace is in the weboptions
                if whitespace: 
                    if line[0:4] == 'ATOM':
                        newline = line[0:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                        pqrfile.write("%s\n" % string.strip(newline))
                    elif line[0:6] == 'HETATM':
                        newline = line[0:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                        pqrfile.write("%s\n" % string.strip(newline))
                else: 
                    pqrfile.write("%s\n" % string.strip(line))
            pqrfile.close()
                    
            if weboptions.otheroptions['apbs']:
                from src import inputgen
                from src import psize
                method = "mg-auto"
                size = psize.Psize()
                size.parseInput(pqrpath)
                size.runPsize(pqrpath)
                async = 0 # No async files here!
                myinput = inputgen.Input(pqrpath, size, method, async)
                myinput.printInputFiles()
                myinput.dumpPickle()
                        
            endtime = time.time() - starttime
            #createResults(header, input, name, endtime, missedligands)
            #logRun(weboptions, endtime, len(lines), weboptions.ff, os.environ["REMOTE_ADDR"])
            #printHeader("PDB2PQR Job Submission",have_opal,jobid=name)
            if 'ligand' in weboptions:
                outputligandfile = open('%s%s%s/%s.mol2' % (INSTALLDIR,TMPDIR, name, name),'w')
                outputligandfile.write(weboptions.ligandfilestring)
                outputligandfile.close()
            outputpdbfile = open('%s%s%s/%s.pdb' % (INSTALLDIR,TMPDIR,name,name),'w')
            outputpdbfile.write(weboptions.pdbfilestring)
            outputpdbfile.close()

            statusfile = open('%s%s%s/pdb2pqr_status' % (INSTALLDIR, TMPDIR, name), 'w')
            statusfile.write('complete\n')
            filelist = glob.glob('%s%s%s/%s*' % (INSTALLDIR, TMPDIR, name, name))
            for filename in filelist:
                statusfile.write(filename+'\n')
            statusfile.close()

    #TODO: Better error reporting.
    #Also, get forked job to properly write error status on failure.
    except StandardError, details:
        print details
        createError(name, details)

def mainCGI():
    """
        Main driver for running PDB2PQR from a web page
    """
    print "Content-type: text/html\n"
    import cgi
    import cgitb

    cgitb.enable()
    form = cgi.FieldStorage()
    
    try:
        weboptions = WebOptions(form)
    except WebOptionsError as error:
        print(error)
        sys.exit(2)
        
    if HAVE_PDB2PQR_OPAL=="1":
        handleOpal(weboptions)
    else:
        handleNonOpal(weboptions)
        
    return
    