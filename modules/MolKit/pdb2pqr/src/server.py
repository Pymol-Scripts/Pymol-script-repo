"""
    CGI Server for PDB2PQR

    This module contains the various functions necessary to run PDB2PQR
    from a web server.

    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Nathan A. Baker (baker@biochem.wustl.edu)
    Todd Dolinsky (todd@ccb.wustl.edu)
    Dept. of Biochemistry and Molecular Biophysics
    Center for Computational Biology
    Washington University in St. Louis

    Jens Nielsen (Jens.Nielsen@ucd.ie)
    University College Dublin

    Additional contributing authors listed in documentation and supporting
    package licenses.

    Copyright (c) 2003-2007.  Washington University in St. Louis.  
    All Rights Reserved.

    This file is part of PDB2PQR.

    PDB2PQR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
 
    PDB2PQR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
   
    You should have received a copy of the GNU General Public License
    along with PDB2PQR; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

    ----------------------------
"""

__date__   = "17 March 2007"
__author__ = "Todd Dolinsky"

import string
import os
import sys
import time

# GLOBAL SERVER VARIABLES

""" The absolute path to root HTML directory """
LOCALPATH   = "/export/home/www/html/pdb2pqr/"

""" The relative path to results directory from script directory.
    The web server (i.e. Apache) MUST be able to write to this directory. """
TMPDIR      = "tmp/"

""" The maximum size of temp directory (in MB) before it is cleaned """
LIMIT       = 500.0

""" The path to the web site *directory* """
WEBSITE     = "http://agave.wustl.edu/pdb2pqr/"

""" The name of the main server page """
WEBNAME     = "server.html"

""" The stylesheet to use """
STYLESHEET  = "http://agave.wustl.edu/css/baker.css"

""" The refresh time (in seconds) for the progress page """
REFRESHTIME = 20

""" The absolute path to the loadavg file - set to "" or None if
    not to be included """
LOADPATH    = "/proc/loadavg"

""" The path to the pdb2pqr log - set to "" or None if not to be
    included.  The web server (i.e. Apache) MUST be able to write to
    this directory. """
LOGPATH     = "%s/%s/usage.txt" % (LOCALPATH, TMPDIR)

def setID(time):
    """
        Given a floating point time.time(), generate an ID.
        Use the tenths of a second to differentiate.

        Parameters
            time:  The current time.time() (float)
        Returns
            id  :  The file id (string)
    """
    strID = "%s" % time
    period = string.find(strID, ".")
    id = "%s%s" % (strID[:period], strID[(period+1):(period+2)])
    return id

def logRun(options, nettime, size, ff, ip):
    """
        Log the CGI run for data analysis.  Log file format is as follows:

        DATE  FF  SIZE  OPTIONS   TIME

        Parameters
            options: The options used for this run (dict)
            nettime: The total time taken for the run (float)
            size:    The final number of non-HETATM atoms in the PDB file (int)
            ff:      The name of the ff used
            ip:      The ip address of the user
    """
    if LOGPATH == "" or LOGPATH == None: return
    date = time.asctime(time.localtime())
    debump  = 0
    opt  = 0
    propka = 0
    file = open(LOGPATH,"a")

    if "ffout" in options: ffout = options["ffout"]
    else: ffout="internal"

    text = "%s\t%s\t%s\t%s\t" % (date, ip, ff, ffout)

    opts = ""
    if "debump" in options:  opts += "debump,"
    if "opt" in options:  opts += "optimize,"
    if "ph" in options: opts += "propka,"
    if "apbs" in options: opts += "apbs,"
    if "chain" in options: opts += "chain,"
    
    if opts == "": opts = "none,"

    text += "%s\t%s\t%.2f\n" % (opts[:-1], size, nettime)
    file.write(text)
    file.close()

def cleanTmpdir():
    """
        Clean up the temp directory for CGI.  If the size of the directory
        is greater than LIMIT, delete the older half of the files.  Since
        the files are stored by system time of creation, this is an
        easier task.
    """
    newdir = []
    size = 0.0
    count = 0
    path = LOCALPATH + TMPDIR

    dir = os.listdir(path)
    for filename in dir:
        size = size + os.path.getsize("%s%s" % (path, filename))
        period = string.find(filename,".")
        id = filename[:period]
        if id not in newdir:
            newdir.append(id)
            count += 1
        
    newdir.sort()
    size = size / (1024.0 * 1024.0)
    
    newcount = 0
    if size >= LIMIT:
        for filename in newdir:
            if newcount > count/2.0: break
            try:
                os.remove("%s%s.pqr" % (path, filename))
            except OSError: pass
            try:
                os.remove("%s%s.in" % (path, filename))
            except OSError: pass
            try:
                os.remove("%s%s.html" % (path, filename))
            except OSError: pass
            newcount += 1

def getQuote(path):
    """
        Get a quote to display for the refresh page.
        Uses fortune to generate a quote.

        Parameters:
            path:   The path to the fortune script (str)
        Returns:
            quote:   The quote to display (str)
    """
    fortune = os.popen(path)
    quote = fortune.read()
    quote = string.replace(quote, "\n", "<BR>")
    quote = string.replace(quote, "\t", "&nbsp;"*5)
    quote = "%s<P>" % quote
    return quote

def printProgress(name, refreshname, reftime, starttime):
    """
        Print the progress of the server

        Parameters
            name:        The ID of the HTML page to write to (string)
            refreshname: The name of the HTML page to refresh to (string)
            reftime:     The length of time to set the refresh wait to (int)
            starttime:   The time as returned by time.time() that the run started (float)
    """
    elapsedtime = time.time() - starttime + REFRESHTIME/2.0 # Add in time offset
    filename = "%s%s%s-tmp.html" % (LOCALPATH, TMPDIR, name)
    file = open(filename,"w")
    file.write("<HTML>\n")
    file.write("<HEAD>\n")
    file.write("<TITLE>PDB2PQR Progress</TITLE>\n")
    file.write("<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET)
    file.write("<meta http-equiv=\"Refresh\" content=\"%s; url=%s\">\n" % \
               (reftime, refreshname))
    file.write("</HEAD>\n")
    file.write("<BODY>\n")
    file.write("<H2>PDB2PQR Progress</H2><P>\n")
    file.write("The PDB2PQR server is generating your results - this page will automatically \n")
    file.write("refresh every %s seconds.<P>\n" % REFRESHTIME)
    file.write("Thank you for your patience!<P>\n")
    file.write("Server Progress:<P>\n")
    file.write("<blockquote>\n")
    file.write("<font size=2>Elapsed Time:</font> <code>%.2f seconds</code><BR>\n" % elapsedtime)
    file.write("</blockquote>\n")
    file.write("Server Information:<P>\n")
    file.write("<blockquote>\n")
    loads = getLoads()
    if loads != None:
        file.write("<font size=2>Server load:</font> <code>%s (1min)  %s (5min)  %s (15min)</code><BR>\n" % (loads[0], loads[1], loads[2]))

    file.write("<font size=2>Server time:</font> <code>%s</code><BR>\n" % (time.asctime(time.localtime())))
    file.write("</blockquote>\n")
    file.write("</BODY></HTML>")
    file.close()

def printAcceptance(name):
    """
        Print the first message to stdout (web browser) - set the
        refresh to the <id>-tmp.html file.

        Parameters
            name:    The ID of the HTML page to redirect to (string)
    """
    waittime = int(REFRESHTIME/2.0)
    print "Content-type: text/html\n"
    print "<HTML>"
    print "<HEAD>"
    print "<TITLE>PDB2PQR Progress</TITLE>"
    print "<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">" % STYLESHEET
    print "<meta http-equiv=\"Refresh\" content=\"%s; url=%s%s%s-tmp.html\">" % \
          (waittime, WEBSITE, TMPDIR, name)
    print "</HEAD>"
    print "<BODY>"
    print "<H2>PDB2PQR Progress</H2><P>"
    print "The PDB2PQR server is generating your results - this page will automatically "
    print "refresh every %s seconds.<P>" % REFRESHTIME
    print "Thank you for your patience!<P>"

    print "Server Information:<P>"
    print "<blockquote>"
    loads = getLoads()
    if loads != None:
        print "<font size=2>Server load:</font> <code>%s (1min)  %s (5min)  %s (15min)</code><BR>" % (loads[0], loads[1], loads[2])

    print "<font size=2>Server time:</font> <code>%s</code><BR>" % (time.asctime(time.localtime()))
    print "</blockquote>"
    
    print "</BODY></HTML>"

def getLoads():
    """
        Get the system load information for output and logging

        Returns
            loads:  A three entry list containing the 1, 5, and
                    15 minute loads. If the load file is not found,
                    return None.
    """
    if LOADPATH == "": return None
    try:
        file = open(LOADPATH)
    except IOError:
        return None

    line = file.readline()
    words = string.split(line)
    loads = words[:3]
    
    return loads

def createResults(header, input, name, time, missedligands=[]):
    """
        Create the results web page for CGI-based runs

        Parameters
            header: The header of the PQR file (string)
            input:   A flag whether an input file has been created (int)
            tmpdir:  The resulting file directory (string)
            name:    The result file root name, based on local time (string)
            time:    The time taken to run the script (float)
            missedligands: A list of ligand names whose parameters could
                     not be assigned. Optional. (list)
    """
    newheader = string.replace(header, "\n", "<BR>")
    newheader = string.replace(newheader," ","&nbsp;")

    filename = "%s%s%s.html" % (LOCALPATH, TMPDIR, name)
    file = open(filename, "w")
    
    file.write("<html>\n")
    file.write("<head>\n")
    file.write("<title>PDB2PQR Results</title>\n")
    file.write("<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET)
    file.write("</head>\n")

    file.write("<body>\n")
    file.write("<h2>PDB2PQR Results</h2>\n")
    file.write("<P>\n")
    file.write("Here are the results from PDB2PQR.  The files will be available on the ")
    file.write("server for a short period of time if you need to re-access the results.<P>\n")
 
    file.write("<a href=\"%s%s%s.pqr\">%s.pqr</a><BR>\n" % (WEBSITE, TMPDIR, name, name))
    if input:
        file.write("<a href=\"%s%s%s.in\">%s.in</a><BR>\n" % (WEBSITE, TMPDIR, name, name))
    pkaname = "%s%s%s.propka" % (LOCALPATH, TMPDIR, name)
    if os.path.isfile(pkaname):
        file.write("<a href=\"%s%s%s.propka\">%s.propka</a><BR>\n" % (WEBSITE, TMPDIR, name, name))
    typename = "%s%s%s-typemap.html" % (LOCALPATH, TMPDIR, name)
    if os.path.isfile(typename):
        file.write("<a href=\"%s%s%s-typemap.html\">%s-typemap.html</a><BR>\n" % (WEBSITE, TMPDIR, name, name)) 
    file.write("<P>The header for your PQR file, including any warnings generated, is:<P>\n")
    file.write("<blockquote><code>\n")
    file.write("%s<P>\n" % newheader)
    file.write("</code></blockquote>\n")
    if missedligands != []:
        file.write("The forcefield that you have selected does not have ")
        file.write("parameters for the following ligands in your PDB file.  Please visit ")
        file.write("<a href=\"http://davapc1.bioch.dundee.ac.uk/programs/prodrg/\">PRODRG</a> ")
        file.write("to convert these ligands into MOL2 format.  This ligand can the be ")
        file.write("parameterized in your PDB2PQR calculation using the PEOE_PB methodology via ")
        file.write("the 'Assign charges to the ligand specified in a MOL2 file' checkbox:<P>\n")
        file.write("<blockquote><code>\n")
        for item in missedligands:
            file.write("%s<BR>\n" % item)
        file.write("<P></code></blockquote>\n")
    file.write("If you would like to run PDB2PQR again, please click <a href=\"%s%s\">\n" % (WEBSITE, WEBNAME))
    file.write("here</a>.<P>\n")
    file.write("<P>Thank you for using the PDB2PQR server!<P>\n")
    file.write("<font size=\"-1\"><P>Total time on server: %.2f seconds</font><P>\n" % time)
    file.write("<font size=\"-1\"><CENTER><I>Last Updated %s</I></CENTER></font>\n" % __date__) 
    file.write("</body>\n")
    file.write("</html>\n")

def createError(name, details):
    """
        Create an error results page for CGI-based runs

        Parameters
            name:    The result file root name, based on local time (string)
            details: The details of the error (string)
    """
    filename = "%s%s%s.html" % (LOCALPATH, TMPDIR, name)
    file = open(filename, "w")

    file.write("<html>\n")
    file.write("<head>\n")
    file.write("<title>PDB2PQR Error</title>\n")
    file.write("<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET)
    file.write("</head>\n")

    file.write("<body>\n")
    file.write("<h2>PDB2PQR Error</h2>\n")
    file.write("<P>\n")
    file.write("An error occurred when attempting to run PDB2PQR:<P>\n")
    file.write("%s<P>\n" % details)
    file.write("If you believe this error is due to a bug, please contact the server administrator.<BR>\n")
    file.write("If you would like to try running PDB2QR again, please click <a href=\"%s%s\">\n" % (WEBSITE, WEBNAME))
    file.write("here</a>.<P>\n")
    file.write("<font size=\"-1\"><CENTER><I>Last Updated %s</I></CENTER></font>\n" % __date__) 
    file.write("</body>\n")
    file.write("</html>\n")
    
def startServer(name):
    """
        Start the PDB2PQR server.  This function is necessary so
        that useful information can be displayed to the user - otherwise
        nothing would be returned until the complete run finishes.

        Parameters
            name:    The ID name of the final file to create (string)
        Returns
            pqrpath: The complete path to the pqr file (string)
    """    
    cleanTmpdir()
    path = LOCALPATH
    tmpdir = TMPDIR

    starttime = time.time()
    pid = os.fork()
    if pid: #Parent - Create refreshed HTML pages and exit
        pid2 = os.fork()
        if pid2: # print to browser and exit
            printAcceptance(name)
            sys.exit()
        else: # Thread in charge of refreshing pages
            home = os.getcwd()
            os.chdir("/")
            os.setsid()
            os.umask(0)
            os.chdir(home)
            os.close(1)
            os.close(2)
            endname = "%s%s%s.html" % (LOCALPATH, TMPDIR, name)
            tmpname = "%s%s%s-tmp.html" % (LOCALPATH, TMPDIR, name)
            while 1:
                if os.path.isfile(endname):
                    refreshname = "%s%s%s.html" % (WEBSITE, TMPDIR, name)
                    reftime = 5
                    printProgress(name, refreshname, reftime, starttime)
                    break
                else:
                    refreshname = "%s%s%s-tmp.html" % (WEBSITE, TMPDIR, name)
                    reftime = REFRESHTIME
                    printProgress(name, refreshname, reftime, starttime)
                time.sleep(REFRESHTIME)

            time.sleep(REFRESHTIME)
            os.remove(tmpname)
            sys.exit()
        
    else: # Child - run PDB2PQR
        home = os.getcwd()
        os.chdir("/")
        os.setsid()
        os.umask(0)
        os.chdir(home)
        os.close(1)
        os.close(2)
        pqrpath = "%s%s%s.pqr" % (path, tmpdir, name)
        return pqrpath
