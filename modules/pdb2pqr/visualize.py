#!@WHICHPYTHON@
"""
    Displays the Jmol input page
"""
print 'Content-type: text/html\n\n'

__date__ = "18 June 2008"
__author__ = "Samir Unni"
__version__ = "0.0.1"

from src.aconf import *
import cgi, cgitb, pickle, urllib, os, glob

cgitb.enable()
form = cgi.FieldStorage()

def initVars():
    #global serviceURL

    if not form.has_key("jobid"):
        pass # add code to redirect to PDB2PQR input page here
    else:
        jobid = form['jobid'].value

        #aspFile = open('./tmp/%s/%s-asp' % (logTime, logTime))
        #appServicePort = pickle.load(aspFile)
        #aspFile.close()

        #jobInfoFile = open('./tmp/%s/%s-jobinfo' % (logTime, logTime))
        #jobID = jobInfoFile.read()
        #jobInfoFile.close()

        aoFile = open('%s%s%s/%s-ao' % (INSTALLDIR, TMPDIR, jobid, jobid))
        apbsOptions = pickle.load(aoFile)
        aoFile.close()

        if APBS_OPAL_URL!="":
            from AppService_client import AppServiceLocator, getOutputsRequest
            apbsOpalJobIDFile = open('%s%s%s/apbs_opal_job_id' % (INSTALLDIR, TMPDIR, jobid))
            apbsOpalJobID = apbsOpalJobIDFile.read()
            apbsOpalJobIDFile.close()
            
            appLocator = AppServiceLocator()
            resp = appLocator.getAppServicePort(APBS_OPAL_URL).getOutputs(getOutputsRequest(apbsOpalJobID))
            if not os.access('%s%s%s' % (INSTALLDIR, TMPDIR, jobid), os.F_OK):
                os.mkdir('%s%s%s' % (INSTALLDIR, TMPDIR, jobid))
            for file in resp._outputFile:
                fileName = file._name
                if fileName!="Standard Output" and fileName!="Standard Error":
                    urllib.urlretrieve(file._url, '%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, fileName))

        return apbsOptions


def main(apbsOptions):
    cgiFile = "jmol.cgi"
    cgiName = "thisform"
    defaultVisType = "jmol"
    checkJmolType = True
    cssFile = 'pdb2pqr.css'
    jobid = form['jobid'].value

    print '<html>'
    print '\t<head>'
    print '\t\t<title>Visualization</title>'
    print '\t\t<link rel="stylesheet" href="pdb2pqr.css" type="text/css">'
    print '\t\t<script type="text/JavaScript" src="jmol/Jmol.js"></script>'
    print '\t\t<script type="text/JavaScript">APPLET_PATH="jmol/";GZIP=""</script>'
    print '\t\t<script type="text/JavaScript" src="jmol/apbsjmol.js"></script>'
    print '\t</head>'
    print '\t<body onload="init()">'
    print '\t\t<script type="text/javascript">createVisualization(%s, -5.0, 5.0)</script>' % (jobid)

    print '\t<script type="text/javascript">'
    print '\tvar gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");'
    print '\tdocument.write(unescape("%3Cscript src=\'" + gaJsHost + "google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E"));'
    print '\t</script>'
    print '\t<script type="text/javascript">'
    print '\ttry {'
    print '\tvar pageTracker = _gat._getTracker("UA-11026338-3");'
    print '\tpageTracker._trackPageview();'
    print '\t} catch(err) {}</script>'
    print '\t</body>'
    print '</html>'


main(initVars())
