#!@WHICHPYTHON@
"""
    APBS input form [fill in later]
"""

__date__ = "22 June 2007"
__author__ = "Samir Unni"
__version__ = "0.0.1"

import string, sys, os, time, errno, shutil, tempfile, urllib, copy, pickle, glob
#sys.path.append("/home/samir/public_html/pdb2pqr")
import src
import cgi, cgitb
import locale
from sys import stdout, stderr, stdin
from src.aconf import *
from src.server import setID
#from initVars import *
#from generateForm import generateForm
#from apbsInputMake import *
#from apbsExec import apbsExec
#from apbsExec import apbsOpalExec
from sgmllib import SGMLParser

def apbsOpalExec(logTime, form, apbsOptions):
    
    sys.path.append(os.path.dirname(HAVE_APBS))
    from ApbsClient import execApbs, initRemoteVars, enoughMemory

    #style = "%spdb2pqr.css" # HARDCODED

    #apbsLog = "apbs_runtime_output.log" 

    #apbsErrs = "apbs_errors.log" 
    
    # Copies PQR file to temporary directory
    pqrFileName = form["pdb2pqrid"].value + '.pqr'
    #shutil.copyfile('../pdb2pqr/tmp/%s' % pqrFileName, './tmp/%s/%s' % (logTime, pqrFileName))
    

    # Removes water from molecule if requested by the user
    if form.has_key("removewater"):
        if form["removewater"].value == "on":
            os.chdir('./tmp/%s' % logTime)
            inpath = pqrFileName 
            infile = open(inpath, "r")
            outpath = inpath[:-4] + '-nowater' + inpath[-4:]
            outfile = open(outpath, "w")
            newinpath = inpath[:-4] + '-water' + inpath[-4:]
            newoutpath = inpath

            while 1:
                line = infile.readline()
                if line == '':
                    break
                if "WAT" in line:
                    pass
                elif "HOH" in line:
                    pass
                else:
                    outfile.write(line)
            infile.close()
            outfile.close()

            shutil.move(inpath, newinpath)
            shutil.move(outpath, newoutpath)
            os.chdir('../../')

    #argv=[os.path.abspath("tmp/%s/%s.in") % (logTime, "apbsinput")] # HARDCODED??  
    argv=[os.path.abspath("%s%s%s/apbsinput.in" % (INSTALLDIR, TMPDIR, logTime))]  
    if DEFAULT_APBS_OPAL_URL == "0":
        vars={'service_url' : APBS_OPAL_URL}
    else:
        vars = None

    # Check for enough memory
    if(not enoughMemory(argv[-1])):
        return 'notenoughmem'

    appServicePortArray = execApbs(vars=vars, argv=argv)

    # if the version number doesn't match, execApbs returns False
    if(appServicePortArray == False):
        return False

    appServicePort = appServicePortArray[0]

    #aspFile = open('./tmp/%s/%s-asp' % (logTime, logTime),'w')
    #pickle.dump(appServicePort, aspFile)
    #aspFile.close()

    resp = appServicePortArray[1]

    return resp._jobID

def apbsExec(logTime, form, apbsOptions):
    
    tempPage = "results.html"

    # Temporary index.py html page - refreshes in 30 seconds
    
    #apbsLog = "apbs_runtime_output.log" 

    #apbsErrs = "apbs_errors.log" 
    
    # Copies PQR file to temporary directory
    pqrFileName = form["pdb2pqrid"].value + '.pqr'
    #shutil.copyfile('../pdb2pqr/tmp/%s' % pqrFileName, './tmp/%s/%s' % (logTime, pqrFileName))
    

    # Removes water from molecule if requested by the user
    try:
        if form["removewater"].value == "on":
            os.chdir('./tmp/%s' % logTime)
            inpath = pqrFileName 
            infile = open(inpath, "r")
            outpath = inpath[:-4] + '-nowater' + inpath[-4:]
            outfile = open(outpath, "w")
            newinpath = inpath[:-4] + '-water' + inpath[-4:]
            newoutpath = inpath

            while 1:
                line = infile.readline()
                if line == '':
                    break
                if "WAT" in line:
                    pass
                elif "HOH" in line:
                    pass
                else:
                    outfile.write(line)
            infile.close()
            outfile.close()

            shutil.move(inpath, newinpath)
            shutil.move(outpath, newoutpath)
            os.chdir('../../')

    except KeyError:
        pass
        
    pid = os.fork()
    if pid:
        print redirector(logTime)
        sys.exit()
    else:
        currentdir = os.getcwd()
        os.chdir("/")
        os.setsid()
        os.umask(0)
        os.chdir(currentdir)
        os.close(1)
        os.close(2)
        #os.chdir('./tmp/%s' % logTime)
        os.chdir('%s%s%s' % (INSTALLDIR, TMPDIR, logTime))
        # LAUNCHING APBS HERE
        statusfile = open('%s%s%s/apbs_status' % (INSTALLDIR, TMPDIR, logTime),'w')
        statusfile.write("running\n")
        statusfile.close()


        apbs_stdin, apbs_stdout, apbs_stderr = os.popen3('%s apbsinput.in' % HAVE_APBS)

        input = open('%s%s%s/apbs_stdout.txt' % (INSTALLDIR, TMPDIR, logTime), 'w')
        input.write(apbs_stdout.read())
        input.close()
        
        endtimefile = open('%s%s%s/apbs_end_time' % (INSTALLDIR, TMPDIR, logTime), 'w')
        endtimefile.write(str(time.time()))
        endtimefile.close()

        input = open('%s%s%s/apbs_stderr.txt' % (INSTALLDIR, TMPDIR, logTime), 'w')
        input.write(apbs_stderr.read())
        input.close()

        statusfile = open('%s%s%s/apbs_status' % (INSTALLDIR, TMPDIR, logTime),'w')
        statusfile.write("complete\n")
        statusfile.write("%s%s%s/apbsinput.in\n" % (INSTALLDIR, TMPDIR, logTime))
        statusfile.write("%s%s%s/%s.pqr\n" % (INSTALLDIR, TMPDIR, logTime, logTime))
        statusfile.write("%s%s%s/io.mc\n" % (INSTALLDIR, TMPDIR, logTime))
        for filename in glob.glob("%s%s%s/%s-*.dx" % (INSTALLDIR, TMPDIR, logTime, logTime)):
            statusfile.write(filename+"\n")
        statusfile.write("%s%s%s/apbs_stdout.txt\n" % (INSTALLDIR, TMPDIR, logTime))
        statusfile.write("%s%s%s/apbs_stderr.txt\n" % (INSTALLDIR, TMPDIR, logTime))
        statusfile.close()
        sys.exit()

def generateForm(file, initVars, pdb2pqrID, type):
    """CGI form generation code"""
    cgifile = "apbs_cgi.cgi"
    cginame = "thisform"
    file.write("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
                <html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
            <head>
                <title>APBS input</title>
                <link href="pdb2pqr.css" type="text/css" rel="stylesheet" />
                <script type=\"text/javascript\">
                    function showHide(elemShow, elemHide) {
                        document.getElementById(elemHide).style.display = 'none';
                        //document.getElementById(elemShow).style.display = 'block';
                        //document.getElementById(elemShow).style.visibility= 'visible';
                        document.getElementById(elemShow).style.display = '';
                    }

                    function toggle(element) {
                        if(document.getElementById(element).style.display == 'none') {
                            document.getElementById(element).style.display = '';
                        }
                        else {
                            document.getElementById(element).style.display = 'none';
                        }
                    }

                    function getElementsByClass(searchClass,node,tag) {
                        var classElements = new Array();
                        if ( node == null )
                            node = document;
                        if ( tag == null )
                            tag = '*';
                        var els = node.getElementsByTagName(tag);
                        var elsLen = els.length;
                        //var pattern = new RegExp("(^|\\\\\\\\s)"+searchClass+"(\\\\\\\\s|$)");
                        for (i = 0, j = 0; i < elsLen; i++) {
                            //if ( pattern.test(els[i].className) ) {
                            if( els[i].className.indexOf(searchClass) != -1) {
                                classElements[j] = els[i];
                                j++;
                            }
                        }
                        return classElements;
                    }

                    function disableCalcType(calcTypeToDisable) {
                        //var elements = getElementsByClass(calcTypeToDisable, null, 'div');
                        var elements = getElementsByClass(calcTypeToDisable);
                        for(var i=0; i<elements.length; i++) {
                            elements[i].style.display = 'none';
                        }
                    }

                    function showCalcType(calcType) {
                        //var calcTypesToDisable = ['mg-auto','mg-para','mg-manual','fe-manual','mg-dummy'];
                        var calcTypesToDisable = ['mg-auto','mg-para','mg-manual','fe-manual','mg-dummy'];
                        for(var i=0; i<calcTypesToDisable.length; i++) {
                            if(calcTypesToDisable[i]!=calcType) {
                                disableCalcType(calcTypesToDisable[i]);
                            }
                        }

                        //var elements = getElementsByClass(calcType, null, 'div');
                        var elements = getElementsByClass(calcType);
                        for(var i=0; i<elements.length; i++) {
                            //elements[i].style.display = 'block';
                            elements[i].style.display = '';
                            //elements[i].style.visibility = 'visible';
                        }
                    }

                    function findCheckedCalc() {
                        for(var i=0; i<document.%s.type.length; i++) {
                            if(document.%s.type[i].checked) {
                                type = document.%s.type[i];
                                return type;
                            }
                        }
                    }

                    function configInitCheckedCalc() {
                        showCalcType('%s');
                    }

                    function executeOnPageLoad(myfunc) {
                        if(window.addEventListener) {
                            window.addEventListener('load', myfunc, false);
                        }

                        else if(window.attachEvent) {
                            window.attachEvent('onload', myfunc);
                        }
                    }

                    executeOnPageLoad(configInitCheckedCalc);


                </script>
            </head>
    <body>
        <!-- ... body of document ... -->
        <h2>APBS web solver</h2>
    """ % (cginame, cginame, cginame, initVars['calculationType'])) # hardcoded css link
    file.write("<h3>Calculation on <a href=\"tmp/%s/%s\" target=\"_blank\">%s</a> with default values provided by PDB2PQR:</h3><br />\n" % (pdb2pqrID, initVars['pqrname'],initVars['pdbID']))
    # Write out the form element
    print "<form action=\"%s\" method=\"post\" enctype=\"multipart/form-data\" name=\"%s\" id=\"%s\">" % (cgifile, cginame, cginame)
    print "<input type=\"submit\" value=\"Launch\"/><br /><br />"
    print """
            If you prefer to run APBS with custom values, click here:
            <input type=\"checkbox\" name=\"customvalues\" onClick=\"toggle(\'params\');"/>
            <br /><br />
            
            <div id=\"params\" style=\"display:none\">
            Please specify the type of calculation (all parameters for the specified calculation must be fullfilled, unless indicated otherwise):
            <br />

    """

    print "<input type=\"radio\" name=\"type\" value=\"mg-auto\" onClick=\"showCalcType(\'mg-auto\');\"",
    if initVars['calculationType'] == "mg-auto":
        print " checked=\"checked\"",
    
    print "/> Automatically-configured sequential focusing multigrid calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#Automatic_sequential_focusing_multigrid_calculation_.28mg-auto.29\" target=\"_blank\"><font title=\"mg-auto\">(<span class=\"tooltip\">?</span>)</font></a>"

    print "<br />"
    print "<input type=\"radio\" name=\"type\" value=\"mg-para\" onClick=\"showCalcType(\'mg-para\');\"",

    if initVars['calculationType'] == "mg-para":
        print " checked=\"checked\"",
    #print "disabled=\"disabled\""

    print """/> Automatically-configured parallel focusing multigrid calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#Automatic_parallel_focusing_multigrid_calculation_.28mg-para.29\" target=\"_blank\"><font title=\"mg-para\">(<span class=\"tooltip\">?</span>)</font></a>"""

    print "<br />"      
    print "<input type=\"radio\" name=\"type\" value=\"mg-manual\" onClick=\"showCalcType(\'mg-manual\');\"",
    if initVars['calculationType'] == "mg-manual":
        print " checked=\"checked\"",
    
    #print "disabled=\"disabled\""

    print """/> Manually-configured multigrid calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#Manual_multigrid_calculation_.28mg-manual.29\" target=\"_blank\"><font title=\"mg-manual\">(<span class=\"tooltip\">?</span>)</font></a>"""


    print "<br />"
    print "<input type=\"radio\" name=\"type\" value=\"fe-manual\" onClick=\"showCalcType(\'fe-manual\');\"",
    if initVars['calculationType'] == "fe-manual":
        print " checked=\"checked\"",
    
    #print "disabled=\"disabled\""

    print """/> Manually-configured adaptive finite element calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#Manual_adaptive_finite_element_calculation_.28fe-manual.29\" target=\"_blank\"><font title=\"fe-manual\">(<span class=\"tooltip\">?</span>)</font></a>"""

    
    print "<br />"

    print "<input type=\"radio\" name=\"type\" value=\"mg-dummy\" onClick=\"showCalcType(\'mg-dummy\');\"",
    if initVars['calculationType'] == "mg-dummy":
        print " checked=\"checked\"",
    
    #print "disabled=\"disabled\""

    print """/> Surface and charge distribution property calculations (does not require solution of the PBE) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#Manual_non-numerical_calculations_.28mg-dummy.29\" target=\"_blank\"><font title=\"mg-dummy\">(<span class=\"tooltip\">?</span>)</font></a>"""

    print "<ul><li><input type=\"checkbox\" name=\"removewater\" value=\"on\" checked=\"checked\"/> Remove water from calculations and visualizations</li></ul>"

    print """
                <div class=\"mg-auto mg-para mg-manual mg-dummy\">
                <ul>"""

    print """
                <table class=\"apbs\" border=\"1\">
                    <tr>
                        <th></th>
                        <th>x-direction</th>
                        <th>y-direction</th>
                        <th>z-direction</th>
                    </tr>
                    <tbody class=\"mg-auto mg-para mg-manual mg-dummy\">
                    <tr>
                        <td>Number of grid points<br />per processor for grid-<br />based discretization <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#dime\" target=\"_blank\"><font title=\"dime\">(<span class=\"tooltip\">?</span>)</font></a></td>
                        <td><input type=\"text\" name=\"dimenx\" size=\"10\" maxlength=\"20\""""
    if initVars.has_key('dime'):
        print "value=\"%d\"" % initVars['dime'][0]
    print "/></td>"
    print "<td><input type=\"text\" name=\"dimeny\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('dime'):
        print "value=\"%d\"" % initVars['dime'][1]
    print "/></td>"
    print "<td><input type=\"text\" name=\"dimenz\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('dime'):
        print "value=\"%d\"" % initVars['dime'][2]
    print "/></td></tr></tbody>"
    # next row
    print "<tbody class=\"mg-auto mg-para mg-dummy\">"
    print "<tr>"
    print "<td>Coarse mesh domain<br />lengths in a focusing<br />calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#cglen\" target=\"_blank\"><font title=\"cglen\">(<span class=\"tooltip\">?</span>)</font></a></td>"
    print "<td><input type=\"text\" name=\"cglenx\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('coarseGridLength'):
        print "value=\"%g\"" % initVars['coarseGridLength'][0]
    print "/>"
    print "<td><input type=\"text\" name=\"cgleny\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('coarseGridLength'):
        print "value=\"%g\"" % initVars['coarseGridLength'][1]
    print "/>"
    print "<td><input type=\"text\" name=\"cglenz\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('coarseGridLength'):
        print "value=\"%g\"" % initVars['coarseGridLength'][2]
    print "/></td></tr></tbody>"
    # next row
    print "<tbody class=\"mg-auto mg-para\">"
    print "<tr>"
    print "<td>Fine mesh domain<br />lengths in a focusing<br />calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#fglen\" target=\"_blank\"><font title=\"fglen\">(<span class=\"tooltip\">?</span>)</font></a></td>"
    print "<td><input type=\"text\" name=\"fglenx\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('fineGridLength'):
        print "value=\"%g\"" % initVars['fineGridLength'][0]
    print "/>"
    print "<td><input type=\"text\" name=\"fgleny\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('fineGridLength'):
        print "value=\"%g\"" % initVars['fineGridLength'][1]
    print "/>"
    print "<td><input type=\"text\" name=\"fglenz\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('fineGridLength'):
        print "value=\"%g\"" % initVars['fineGridLength'][2]
    print "/></td></tr></tbody>"

    # next row
    print "<tbody class=\"mg-para\">"
    print "<tr>"
    print "<td>Number of proces-<br />sors in a parallel<br />focusing calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#pdime\" target=\"_blank\"><font title=\"pdime\">(<span class=\"tooltip\">?</span>)</font></a></td>"
    print "<td><input type=\"text\" name=\"pdimex\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('pdime'):
        print "value=\"%g\"" % initVars['pdime'][0]
    print "/>"
    print "<td><input type=\"text\" name=\"pdimey\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('pdime'):
        print "value=\"%g\"" % initVars['pdime'][1]
    print "/>"
    print "<td><input type=\"text\" name=\"pdimez\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('pdime'):
        print "value=\"%g\"" % initVars['pdime'][2]
    print "/></td></tr></tbody>"

    # next row
    print "<tbody class=\"mg-manual\">"
    print "<tr>"
    print "<td>Mesh domain<br />lengths <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#glen\" target=\"_blank\"><font title=\"glen\">(<span class=\"tooltip\">?</span>)</font></a></td>"
    print "<td><input type=\"text\" name=\"glenx\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('glen'):
        print "value=\"%g\"" % initVars['glen'][0]
    print "/>"
    print "<td><input type=\"text\" name=\"gleny\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('glen'):
        print "value=\"%g\"" % initVars['glen'][1]
    print "/>"
    print "<td><input type=\"text\" name=\"glenz\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('glen'):
        print "value=\"%g\"" % initVars['glen'][2]
    print "/></td></tr></tbody>"



    print "</table></div></ul>"

    print """           <div class=\"mg-para\"><ul>
                <li>Amount of overlap to include between the individual processors' meshes <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ofrac\" target=\"_blank\"><font title=\"frac\">(<span class=\"tooltip\">?</span>)</font></a>:"""
    print "<input type=\"text\" name=\"frac\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('processorMeshOverlap'):
        print "value=\"%f\"" % initVars['processorMeshOverlap']
    print "/></li><br />"

    print "<li><input type=\"checkbox\" name=\"asyncflag\" onClick=toggle(\"async\") "
    if initVars['asyncflag']:
        print " checked=\"checked\""
    print """/> Perform the tasks in a parallel run asynchronously <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#async\" target=\"_blank\"><font title=\"asyncflag\">(<span class=\"tooltip\">?</span>)</font></a></li>"""
    print "<blockquote>"
    print "<div id=\"async\" style=\"display:none\">"
    print "<li>Rank for a processor to masquerade as <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#async\" target=\"_blank\"><font title=\"async\">(<span class=\"tooltip\">?</span>)</font></a>:"
    print "<input type=\"text\" name=\"async\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('async'):
        print " value=\"%i\"" % initVars['async']
    print "/></li>"
    print "</blockquote>"


    print "</li></ul></div>"

    print "<div class=\"mg-manual mg-dummy\">"
    print "<ul><li>Depth of the multilevel hierarchy used in the multigrid solver <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#nlev\" target=\"_blank\"><font title=\"nlev\">(<span class=\"tooltip\">?</span>)</font></a>:"
    print "<input type=\"text\" name=\"nlev\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('nlev'):
        print " value=\"%i\"" % initVars['nlev']
    print "/></li></ul></div>"

    print """           <div class=\"mg-manual mg-dummy\"><ul>
                
                <li>Center of the grid <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#gcent\" target=\"_blank\"><font title=\"gcent\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>
                <blockquote><ul><li>"""
    print "<input type=\"radio\" name=\"gcent\" value=\"mol\" onClick=\"showHide(\'gcentmol\',\'gcentcoord\');\""

    if initVars['gridCenterMethod'] == "molecule":
        print "checked=\"checked\""

    print "/> Center the grid on a molecule.</li></ul>"
    
    print "<div id=\"gcentmol\""
    if initVars['gridCenterMethod'] != "molecule":
        print " style=\"display: none;\""
    print"""><blockquote>
                    <ul>
                        <li>Enter molecule ID:"""
    
    print "<input type=\"text\" name=\"gcentid\" size=\"10\" maxlength=\"20\""
    if initVars['gridCenterMethod'] == "molecule" and initVars.has_key('gridCenterMoleculeID'):
        print "value=\"%d\"" % initVars['gridCenterMoleculeID']
    print "/>"

    print """
                        </li>
                    </ul>
                    </blockquote></div>
                    """

    print "<ul><li><input type=\"radio\" name=\"gcent\" value=\"coord\" onClick=\"showHide(\'gcentcoord\',\'gcentmol\');\""

    if initVars['gridCenterMethod'] == "coordinate":
        print "checked"

    print "/> Manually enter coordinates for center of grid:</li></ul>"

    print "<div id=\"gcentcoord\""
    if initVars['gridCenterMethod'] != "coordinate":
        print "style=\"display: none;\""
    print """><blockquote>
                    <ul>
                        <li>"""

    print "x-coordinate: <input type=\"text\" name=\"gxcent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('gridCenter'):
        print "value=\"%d\"" % initVars['gridCenter'][0]
    print "/>"


    print """
                        </li>
                        <li>"""

    print "y-coordinate: <input type=\"text\" name=\"gycent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('gridCenter'):
        print "value=\"%d\"" % initVars['gridCenter'][1]
    print "/>"

    print """
                        </li>
                        <li>"""

    print "z-coordinate: <input type=\"text\" name=\"gzcent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('gridCenter'):
        print "value=\"%d\"" % initVars['gridCenter'][2]
    print "/>"
    print "</li></ul></blockquote></div>"
    print "</blockquote>"
    print "</div>"












    print """ 
        
       <div class=\"mg-auto mg-para\"><ul>
                <li>Center of the coarse grid <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#cgcent\" target=\"_blank\"><font title=\"cgcent\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>
                                <blockquote><ul><li>"""

    print "<input type=\"radio\" name=\"cgcent\" value=\"mol\" onClick=\"showHide(\'cgcentmol\',\'cgcentcoord\');\""

    if initVars['coarseGridCenterMethod'] == "molecule":
        print "checked=\"checked\""

    print "/> Center the grid on a molecule.</li></ul>"
    
    print "<div id=\"cgcentmol\""
    if initVars['coarseGridCenterMethod'] != "molecule":
        print " style=\"display: none;\""
    print"""><blockquote>
                    <ul>
                        <li>Enter molecule ID:"""
    
    print "<input type=\"text\" name=\"cgcentid\" size=\"10\" maxlength=\"20\""
    if initVars['coarseGridCenterMethod'] == "molecule" and initVars.has_key('coarseGridCenterMoleculeID'):
        print "value=\"%d\"" % initVars['coarseGridCenterMoleculeID']
    print "/>"

    print """
                        </li>
                    </ul>
                    </blockquote>
                    </div>
                    """

    print "<ul><li><input type=\"radio\" name=\"cgcent\" value=\"coord\" onClick=\"showHide(\'cgcentcoord\',\'cgcentmol\');\""

    if initVars['coarseGridCenterMethod'] == "coordinate":
        print "checked"

    print "/> Manually enter coordinates for center of grid:</li></ul>"

    print "<div id=\"cgcentcoord\""
    if initVars['coarseGridCenterMethod'] != "coordinate":
        print "style=\"display: none;\""
    print """><blockquote>
                    <ul>
                        <li>"""

    print "x-coordinate: <input type=\"text\" name=\"cgxcent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('coarseGridCenter'):
        print "value=\"%d\"" % initVars['coarseGridCenter'][0]
    print "/>"


    print """
                        </li>
                        <li>"""

    print "y-coordinate: <input type=\"text\" name=\"cgycent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('coarseGridCenter'):
        print "value=\"%d\"" % initVars['coarseGridCenter'][1]
    print "/>"

    print """
                        </li>
                        <li>"""

    print "z-coordinate: <input type=\"text\" name=\"cgzcent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('coarseGridCenter'):
        print "value=\"%d\"" % initVars['coarseGridCenter'][2]
    print "/>"
    print """</li></ul></blockquote>
                    </blockquote>"""
            #"""</div>"""

    print """


                <ul>
                <li>Center of the fine grid <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#fgcent\" target=\"_blank\"><font title=\"fgcent\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>
                <blockquote>
                <ul><li>"""

    print "<input type=\"radio\" name=\"fgcent\" value=\"mol\" onClick=\"showHide(\'fgcentmol\',\'fgcentcoord\');\""

    if initVars['fineGridCenterMethod'] == "molecule":
        print "checked=\"checked\""

    print "/> Center the grid on a molecule.</li></ul>"
    print "<div id=\"fgcentmol\""
    if initVars['fineGridCenterMethod'] != "molecule":
        print "style=\"display: none\""

    print """
                    ><blockquote>
                    <ul>
                        <li>Enter molecule ID:"""
    
    print "<input type=\"text\" name=\"fgcentid\" size=\"10\" maxlength=\"20\""
    if initVars['fineGridCenterMethod'] == "molecule" and initVars.has_key('fineGridCenterMoleculeID'):
        print "value=\"%d\"" % initVars['fineGridCenterMoleculeID']
    print "/>"
    
    print """
                        </li>
                    </ul>
                    </blockquote></div>"""
    
    print "<ul><li><input type=\"radio\" name=\"fgcent\" value=\"coord\" onClick=\"showHide(\'fgcentcoord\',\'fgcentmol\');\""
    if initVars['fineGridCenterMethod'] == "coordinate":
        print "checked"
    print "/> Manually enter coordinates for the center of the grid.</li></ul>"
    print "<div id=\"fgcentcoord\""
    if initVars['fineGridCenterMethod'] != "coordinate":
        print "style=\"display: none;\""

    print """
                    ><blockquote>
                    <ul>
                        <li>"""

    print "x-coordinate: <input type=\"text\" name=\"fgxcent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('fineGridCenter'):
        print "value=\"%d\"" % initVars['fineGridCenter'][0]
    print "/>"
    print """
                        </li>
                        <li>"""

    print "y-coordinate: <input type=\"text\" name=\"fgycent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('fineGridCenter'):
        print "value=\"%d\"" % initVars['fineGridCenter'][1]
    print "/>"

    print """
                        </li>
                        <li>"""

    print "z-coordinate: <input type=\"text\" name=\"fgzcent\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('fineGridCenter'):
        print "value=\"%d\"" % initVars['fineGridCenter'][2]
    print "/>"
    
    print """
                        </li>
                </ul>
                    </blockquote>
                </blockquote>
                    </div>"""#</div>"""


    #print       """<ul>
    #            <div class=\"mg-para mg-manual fe-manual mg-dummy\""""
    #print """>
    #            <li>Molecule for which the PBE is to be solved <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#mol\" target=\"_blank\"><font title=\"mol\">(<span class=\"tooltip\">?</span>)</font></a>:"""

    #print "<input type=\"text\" name=\"mol\" size=\"10\" maxlength=\"20\""
    #if initVars.molecule != None:
    #    print "value=\"%d\"" % initVars.molecule
    #print "/></li></div>"

    #print """</ul>


    print """   <ul>
                <li>Type of PBE to be solved:</li></ul>"""

    print """<blockquote>
    <ul>"""
    print "<li><input type=\"radio\" name=\"solvetype\" value=\"lpbe\""
    if initVars['solveType'] == "linearized":
        print "checked=\"checked\""
    print "/> Linearized <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#lpbe\" target=\"_blank\"><font title=\"lpbe\">(<span class=\"tooltip\">?</span>)</font></a></li>"


    print "<li><input type=\"radio\" name=\"solvetype\" value=\"npbe\""
    if initVars['solveType'] == "nonlinearized":
        print "checked=\"checked\""
    print "/> Nonlinearized <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#npbe\" target=\"_blank\"><font title=\"npbe\">(<span class=\"tooltip\">?</span>)</font></a></li>"

    print "<div class=\"fe-manual\""
    #if initVars.defaultCalcType != "fe-manual":
    #    print " style=\"display: none;\""
    print "><li><input type=\"radio\" name=\"solvetype\" value=\"lrpbe\""
    if initVars['solveType'] == "linearized regularized":
        print "checked=\"checked\""
    print "/> Linearized (regularized) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#lrpbe\" target=\"_blank\"><font title=\"lrpbe\">(<span class=\"tooltip\">?</span>)</font></a></li>"

    
    print "<li><input type=\"radio\" name=\"solvetype\" value=\"nrpbe\""
    if initVars['solveType'] == "nonlinearized regularized":
        print "checked=\"checked\""
    print "/> Nonlinearized (regularized) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#nrpbe\" target=\"_blank\"><font title=\"nrpbe\">(<span class=\"tooltip\">?</span>)</font></a></li></div>"

    print """</ul>
    </blockquote>"""

    print "         <ul><li>Boundary condition definition <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#bcfl\" target=\"_blank\"><font title=\"bcfl\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"
    
    print "<blockquote> <ul>"
    print "<li><input type=\"radio\" name=\"bcfl\" value=\"zero\""
    if initVars['boundaryConditions'] == "zero":
        print "checked=\"checked\""
    print "/> Zero <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#bcfl\" target=\"_blank\"><font title=\"zero\">(<span class=\"tooltip\">?</span>)</font></a></li>"


    print "<li><input type=\"radio\" name=\"bcfl\" value=\"sdh\""
    if initVars['boundaryConditions'] == "sdh":
        print "checked=\"checked\""
    print "/> Single Debye-Huckel <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#bcfl\" target=\"_blank\"><font title=\"sdh\">(<span class=\"tooltip\">?</span>)</font></a></li>"


    print "<li><input type=\"radio\" name=\"bcfl\" value=\"mdh\""
    if initVars['boundaryConditions'] == "mdh":
        print "checked=\"checked\""
    print "/> Multiple Debye-Huckel <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#bcfl\" target=\"_blank\"><font title=\"mdh\">(<span class=\"tooltip\">?</span>)</font></a></li>"


    print "<li><input type=\"radio\" name=\"bcfl\" value=\"focus\""
    if initVars['boundaryConditions'] == "focus":
        print "checked=\"checked\""
    print "/> Focusing <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#bcfl\" target=\"_blank\"><font title=\"focus\">(<span class=\"tooltip\">?</span>)</font></a></li>"
    print "</ul></blockquote>"


    print "<ul>"
    print """       <li>Mobile ion species present in system (optional) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"ion\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""

    print """<ul>
                <table class=\"apbs\" border=\"1\">
                    <tr>
                        <th></th>
                        <th>Mobile ion species<br />charge (in e<sub>c</sub>) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"charge\">(<span class=\"tooltip\">?</span>)</font></a></th>
                        <th>Mobile ion species<br />concentration (in M) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"conc\">(<span class=\"tooltip\">?</span>)</font></a></th>
                        <th>Mobile ion species<br />radius (in A) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"radius\">(<span class=\"tooltip\">?</span>)</font></a></th>
                    </tr>"""
    # new row
    for i in range(0,3):
        print """<tr>
                            <td>Ion %d</td>
                            <td><input type=\"text\" name=\"charge%d\" size=\"10\" maxlength=\"20\"""" % ((i+1),i)
        if initVars.has_key('mobileIonSpeciesCharge'):
            print "value=\"%d\"" % initVars['mobileIonSpeciesCharge']
        print "/></td>"
        print """
                            <td><input type=\"text\" name=\"conc%d\" size=\"10\" maxlength=\"20\"""" % i
        if initVars.has_key('mobileIonSpeciesConcentration'):
            print "value=\"%d\"" % initVars['mobileIonSpeciesConcentration']
        print "/></td>"
        print """
                            <td><input type=\"text\" name=\"radius%d\" size=\"10\" maxlength=\"20\"""" % i
        if initVars.has_key('mobileIonSpeciesRadius'):
            print "value=\"%d\"" % initVars['mobileIonSpeciesRadius']
        print "/></td>"
    print "</tr></table></ul>"

                        

            
    #print """<blockquote><ul>
    #            <li>Mobile ion species charge (in e<sub>c</sub>) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"charge\">(<span class=\"tooltip\">?</span>)</font></a>: """

    #print "<input type=\"text\" name=\"charge\" size=\"10\" maxlength=\"20\""
    #if initVars.defaultMobileIonSpeciesCharge != None:
    #    print "value=\"%d\"" % initVars.defaultMobileIonSpeciesCharge
    #print "/></li>"

    #
    #print "         <li>Mobile ion species concentration (in M) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"conc\">(<span class=\"tooltip\">?</span>)</font></a>:"
    #
    #print "<input type=\"text\" name=\"conc\" size=\"10\" maxlength=\"20\""
    #if initVars.defaultMobileIonSpeciesConcentration != None:
    #    print "value=\"%d\"" % initVars.defaultMobileIonSpeciesConcentration
    #print "/></li>"

    #print "         <li>Mobile ion species radius (in A) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#ion\" target=\"_blank\"><font title=\"radius\">(<span class=\"tooltip\">?</span>)</font></a>:"
    #
    #print "<input type=\"text\" name=\"radius\" size=\"10\" maxlength=\"20\""
    #if initVars.defaultMobileIonSpeciesRadius != None:
    #    print "value=\"%d\"" % initVars.defaultMobileIonSpeciesRadius
    #print "/></li>"
    #print "</ul></blockquote>"


    print "         <ul><li>Biomolecular dielectric constant <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#pdie\" target=\"_blank\"><font title=\"pdie\">(<span class=\"tooltip\">?</span>)</font></a>:"

    
    print "<input type=\"text\" name=\"pdie\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('biomolecularDielectricConstant'):
        print "value=\"%g\"" % initVars['biomolecularDielectricConstant']
    print "/>"
    
    print """       </li>
                </ul>"""

    print "         <ul><li>Dielectric constant of solvent <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#sdie\" target=\"_blank\"><font title=\"sdie\">(<span class=\"tooltip\">?</span>)</font></a>:"

    
    print "         <input type=\"text\" name=\"sdie\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('dielectricSolventConstant'):
        print "value=\"%g\"" % initVars['dielectricSolventConstant']
    
    print """/></li>
                </ul>
    
                <ul>
                <li>Method by which the biomolecular point charges are mapped onto the grid <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#chgm\" target=\"_blank\"><font title=\"chgm\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""

    print "<blockquote><ul>"
    print "<li><input type=\"radio\" name=\"chgm\" value=\"spl0\""
    if initVars.has_key('biomolecularPointChargeMapMethod') and initVars['biomolecularPointChargeMapMethod'] == "spl0":
        print "checked"
    print "/> Traditional trilinear interpolation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#chgm\" target=\"_blank\"><font title=\"spl0\">(<span class=\"tooltip\">?</span>)</font></a></li>"
    

    print "<li><input type=\"radio\" name=\"chgm\" value=\"spl2\""
    if initVars.has_key('biomolecularPointChargeMapMethod') and initVars['biomolecularPointChargeMapMethod'] == "spl2":
        print "checked=\"checked\""
    print "/> Cubic B-spline discretization <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#chgm\" target=\"_blank\"><font title=\"spl2\">(<span class=\"tooltip\">?</span>)</font></a></li>"

    
    print "<li><input type=\"radio\" name=\"chgm\" value=\"spl4\""
    if initVars.has_key('biomolecularPointChargeMapMethod') and initVars['biomolecularPointChargeMapMethod'] == "spl4":
        print "checked"
    print "/> Quintic B-spline discretization <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#chgm\" target=\"_blank\"><font title=\"spl4\">(<span class=\"tooltip\">?</span>)</a></font></li>"
    print "</ul></blockquote>"


    print "         <ul><li>Number of grid points per square-angstrom to use in surface constructions <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#sdens\" target=\"_blank\"><font title=\"sdens\">(<span class=\"tooltip\">?</span>)</font></a>:"
    
    print "<input type=\"text\" name=\"sdens\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('surfaceConstructionResolution'):
        print "value=\"%g\"" % initVars['surfaceConstructionResolution']

    print """   /></li></ul>


            <ul><li>Model to use to construct the dielectric ion-accessibility coefficients <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#srfm\" target=\"_blank\"><font title=\"srfm\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""

    
    print "<blockquote> <ul>"
    print "<li><input type=\"radio\" name=\"srfm\" value=\"mol\""
    if initVars.has_key('dielectricIonAccessibilityModel') and initVars['dielectricIonAccessibilityModel'] == "mol":
        print "checked=\"checked\""
    
    print "/> Dielectric coefficient is defined based on a molecular surface definition; ion-accessibility coefficient is defined by an \"inflated\" van der Waals model <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#srfm\" target=\"_blank\"><font title=\"mol\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"

    print "<li><input type=\"radio\" name=\"srfm\" value=\"smol\""
    if initVars.has_key('dielectricIonAccessibilityModel') and initVars['dielectricIonAccessibilityModel'] == "smol":
        print "checked=\"checked\""
    
    print "/> Dielectric and ion-accessiblity coefficients are defined as above, but then are then \"smoothed\" by a 9-point harmonic averaging to somewhat reduce sensitivity to the grid setup <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#srfm\" target=\"_blank\"><font title=\"smol\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"

    print "<li><input type=\"radio\" name=\"srfm\" value=\"spl2\""
    if initVars.has_key('dielectricIonAccessibilityModel') and initVars['dielectricIonAccessibilityModel'] == "spl2":
        print "checked=\"checked\""
    
    print "/> Dielectric and ion-accessibility coefficients are defined by a cubic-spline surface <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#srfm\" target=\"_blank\"><font title=\"spl2\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"

    print "<li><input type=\"radio\" name=\"srfm\" value=\"spl4\""
    if initVars.has_key('dielectricIonAccessibilityModel') and initVars['dielectricIonAccessibilityModel'] == "spl4":
        print "checked=\"checked\""

    print "/> Dielectric and ion-accessibility coefficients are defined by a 7th order polynomial <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#srfm\" target=\"_blank\"><font title=\"spl4\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"
    print "</ul></blockquote>"


    print "         <ul><li>Radius of the solvent molecules <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#srad\" target=\"_blank\"><font title=\"srad\">(<span class=\"tooltip\">?</span>)</font></a>:"

    
    print "<input type=\"text\" name=\"srad\" size=\"10\" maxlength=\"20\""
    if initVars['solventRadius'] != None:
        print "value=\"%g\"" % initVars['solventRadius']
    
    print """/></li></ul>


                <ul><li>Size of the support for spline-based surface definitions <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#swin\" target=\"_blank\"><font title=\"swin\">(<span class=\"tooltip\">?</span>)</font></a>:"""
    
    
    print "<input type=\"text\" name=\"swin\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('surfaceDefSupportSize'):
        print "value=\"%g\"" % initVars['surfaceDefSupportSize']
    
    print """/></li></ul>


                <ul><li>Temperature for PBE calculation (in K) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#temp\" target=\"_blank\"><font title=\"temp\">(<span class=\"tooltip\">?</span>)</font></a>:"""
    
    
    print "<input type=\"text\" name=\"temp\" size=\"10\" maxlength=\"20\""
    if initVars.has_key('temperature'):
        print "value=\"%g\"" % initVars['temperature']
    
    print """/></li></ul>


                <ul><li>Calculation of electrostatic energy from a PBE calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcenergy\" target=\"_blank\"><font title=\"calcenergy\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""
    
    
    print "<blockquote><ul>"
    print "<li><input type=\"radio\" name=\"calcenergy\" value=\"no\""
    if initVars['calculationEnergy'] == "no":
        print "checked"
    
    print """/> Don\'t calculate any energies <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcenergy\" target=\"_blank\"><font title=\"no\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    print "<li><input type=\"radio\" name=\"calcenergy\" value=\"total\""
    if initVars['calculationEnergy'] == "total":
        print "checked=\"checked\""
    
    print """/> Calculate and return total electrostatic energy for the entire molecule <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcenergy\" target=\"_blank\"><font title=\"total\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    print "<li><input type=\"radio\" name=\"calcenergy\" value=\"comps\""
    if initVars['calculationEnergy'] == "comps":
        print "checked=\"checked\""
    
    print """/> Calculate and return total electrostatic energy for the entire molecule as well as electrostatic energy components for each atom <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcenergy\" target=\"_blank\"><font title=\"comps\">(<span class=\"tooltip\">?</span>)</font></a>.</li>
    
    </ul></blockquote>


                <ul><li>Calculation of electrostatic and apolar force outputs from a PBE calculation <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcforce\" target=\"_blank\"><font title=\"calcforce\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""

    
    print "<blockquote><ul>"
    print "<li><input type=\"radio\" name=\"calcforce\" value=\"no\""
    if initVars['calculationForce'] == "no":
        print "checked=\"checked\""
    
    print """/> Don\'t calculate any forces <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcforce\" target=\"_blank\"><font title=\"no\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    print "<li><input type=\"radio\" name=\"calcforce\" value=\"total\""
    if initVars['calculationForce'] == "total":
        print "checked=\"checked\""
    
    print """/> Calculate and return total electrostatic and apolar forces for the entire molecule <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcforce\" target=\"_blank\"><font title=\"total\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    print "<li><input type=\"radio\" name=\"calcforce\" value=\"comps\""
    if initVars['calculationForce'] == "comps":
        print "checked=\"checked\""
    
    print """/> Calculate and return total electrostatic and apolar forces for the entire molecule as well as force components for each atom <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#calcforce\" target=\"_blank\"><font title=\"comps\">(<span class=\"tooltip\">?</span>)</font></a>.</li>
    </ul> </blockquote>"""

    print """       

                <ul><li>Output of scalar data calculated during the PB run <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"write\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""

    
    print "<blockquote><ul>"
    print "<li><input type=\"checkbox\" name=\"writecharge\""
    if initVars['writeBiomolecularChargeDistribution']:
        print "checked=\"checked\""

    print """/> Write out the biomolecular charge distribution in units of e<sub>c</sub> (multigrid only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"charge\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writepot\""
    if initVars['writeElectrostaticPotential']:
        print "checked=\"checked\""

    print """/> Write out the electrostatic potential in units of k<sub>b</sub>T/e<sub>c</sub>  (multigrid and finite element) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"pot\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writesmol\""
    if initVars['writeMolecularSurfaceSolventAccessibility']:
        print "checked=\"checked\""
    
    print """/> Write out the solvent accessibility defined by the molecular surface definition <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"smol\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    
    print "<li><input type=\"checkbox\" name=\"writesspl\""
    if initVars['writeSplineBasedSolventAccessibility']:
        print "checked=\"checked\""

    print """/> Write out the spline-based solvent accessibility <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"sspl\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    
    print "<li><input type=\"checkbox\" name=\"writevdw\""
    if initVars['writeVanDerWaalsSolventAccessibility']:
        print "checked=\"checked\""
    
    print """/> Write out the van der Waals-based solvent accessibility <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"blank\"><font title=\"vdw\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writeivdw\""
    if initVars['writeInflatedVanDerWaalsIonAccessibility']:
        print "checked=\"checked\""
    
    print """/> Write out the inflated van der Waals-based ion accessibility <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"ivdw\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""
                
    
    print "<li><input type=\"checkbox\" name=\"writelap\""
    if initVars['writePotentialLaplacian']:
        print "checked=\"checked\""
    
    print """/> Write out the Laplacian of the potential in units of k<sub>B</sub>T/e<sub>c</sub>/A<sup>2</sup> (multigrid only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"lap\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writeedens\""
    if initVars['writeEnergyDensity']:
        print "checked=\"checked\""
    
    print """/> Write out the \"energy density\" in units of k<sub>B</sub>T/e<sub>c</sub>/A<sup>2</sup> (multigrid only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"edens\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

    
    print "<li><input type=\"checkbox\" name=\"writendens\""
    if initVars['writeMobileIonNumberDensity']:
        print "checked=\"checked\""
    
    print """/> Write out the mobile ion number density for <i>m</i> ion species in units of M (multigrid only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"ndens\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writeqdens\""
    if initVars['writeMobileChargeDensity']:
        print "checked=\"checked\""
    
    print """/> Write out the mobile charge density for <i>m</i> ion species in units of e<sub>c</sub> M (multigrid only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"qdens\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writedielx\""
    if initVars['writeDielectricMapShift'][0]:
        print "checked=\"checked\""
    
    print """/> Write out the dielectric map shifted by <sup>1</sup>/<sub>2</sub> grid spacing in the x-direction <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"dielx\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writediely\""
    if initVars['writeDielectricMapShift'][1]:
        print "checked=\"checked\""
    
    print """/> Write out the dielectric map shifted by <sup>1</sup>/<sub>2</sub> grid spacing in the y-direction <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"diely\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writedielz\""
    if initVars['writeDielectricMapShift'][2]:
        print "checked=\"checked\""
    
    print """/> Write out the dielectric map shifted by <sup>1</sup>/<sub>2</sub> grid spacing in the z-direction <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"dielz\">(<span class=\"tooltip\">?</span>)</font></a>.</li>"""

                
    print "<li><input type=\"checkbox\" name=\"writekappa\""
    if initVars['writeIonAccessibilityKappaMap']:
        print "checked=\"checked\""
    
    print """/> Write out the ion-accessibility kappa map <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"kappa\">(<span class=\"tooltip\">?</span>)</font></a>.</li>
    </ul></blockquote>


                <ul><li>Format for writing out the data <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"format\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""

    print "<blockquote><ul>"
    print "<li><input type=\"radio\" name=\"writeformat\" value=\"dx\""
    if initVars['format'] == "dx":
        print "checked=\"checked\""
    
    print """/> OpenDX (multigrid and finite element) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"dx\">(<span class=\"tooltip\">?</span>)</font></a></li>"""
    
    print """
                <li><input type=\"radio\" name=\"writeformat\" value=\"avs\" disabled=\"disabled\"/> AVS UCD (finite element only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"avs\">(<span class=\"tooltip\">?</span>)</font></a></li>
                <li><input type=\"radio\" name=\"writeformat\" value=\"uhbd\" disabled=\"disabled\"/> UBHD (multigrid only) <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#write\" target=\"_blank\"><font title=\"uhbd\">(<span class=\"tooltip\">?</span>)</font></a></li>
                </ul></blockquote>
            </div>"""



    #print """<ul><li>Choose type of operator to output <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#writemat\" target=\"_blank\"><font title=\"type\">(<span class=\"tooltip\">?</span>)</font></a>:</li></ul>"""
                 #ADD PYTHON CODE ABOVE WHEN OTHER OPTIONS ARE AVAILABLE
                
    #print "<blockquote><ul>"
    #print """           <li><input type=\"radio\" name=\"writemat\" value=\"poisson\" disabled=\"disabled\"/> Poisson <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#writemat\" target=\"_blank\"><font title=\"poisson\">(<span class=\"tooltip\">?</span>)</font></a></li>
                #<li><input type=\"radio\" name=\"writemat\" value=\"pot\" disabled=\"disabled\"/> Gateaux derivative of the full PBE operator evaluated at the current solution <a href=\"http://apbs.wustl.edu/MediaWiki/index.php/ELEC_input_file_section#writemat\" target=\"_blank\"><font title=\"pot\">(<span class=\"tooltip\">?</span>)</font></a></li>
                #</ul>
                #</blockquote>
                
                #""" #ADD PYTHON CODE ABOVE WHEN OPTIONS ARE AVAILABLE

    print "<br />"
    if type=="local":
        print "<input type=\"hidden\" name=\"hiddencheck\" value=\"local\"/>"
    else:
        print "<input type=\"hidden\" name=\"hiddencheck\" value=\"opal\"/>"

    print "<input type=\"hidden\" name=\"pdb2pqrid\"",
    print "value=\"%s\"" % pdb2pqrID,
    print "/>"

    print "<input type=\"hidden\" name=\"mol\" value=\"1\"/>"

    print """
        </form> 
    <p>
        <a href="http://validator.w3.org/check?uri=referer"><img
                src="http://www.w3.org/Icons/valid-xhtml10"
                    alt="Valid XHTML 1.0 Transitional" height="31" width="88" /></a>
                  </p>
<script type="text/javascript">
var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
try {
var pageTracker = _gat._getTracker("UA-11026338-3");
pageTracker._trackPageview();
} catch(err) {}</script>
    </body>
</html>

        """     

def unpickleVars(pdb2pqrID):
    """ Converts instance pickle from PDB2PQR into a dictionary for APBS """
    apbsOptions = {}
    #pfile = open("/home/samir/public_html/pdb2pqr/tmp/%s-input.p" % pdb2pqrID, 'r')
    pfile = open("%s%s%s/%s-input.p" % (INSTALLDIR, TMPDIR, pdb2pqrID, pdb2pqrID), 'r')
    inputObj = pickle.load(pfile)
    pfile.close()
    myElec = inputObj.elecs[0]

    apbsOptions['pqrname'] = pdb2pqrID+'.pqr'
    apbsOptions['pdbID'] = inputObj.pqrname[:-4]
    
    if myElec.cgcent[0:3] == "mol":
        apbsOptions['coarseGridCenterMethod'] = "molecule"
        apbsOptions['coarseGridCenterMoleculeID'] = locale.atoi(myElec.cgcent[4:])
    else:
        apbsOptions['coarseGridCenterMethod'] = "coordinate"
        apbsOptions['coarseGridCenter'] = myElec.cgcent

    if myElec.fgcent[0:3] == "mol":
        apbsOptions['fineGridCenterMethod'] = "molecule"
        apbsOptions['fineGridCenterMoleculeID'] = locale.atoi(myElec.fgcent[4:])
    else:
        apbsOptions['fineGridCenterMethod'] = "coordinate"
        apbsOptions['fineGridCenter'] = myElec.fgcent

    if myElec.gcent[0:3] == "mol":
        apbsOptions['gridCenterMethod'] = "molecule"
        apbsOptions['gridCenterMoleculeID'] = locale.atoi(myElec.gcent[4:])
    else:
        apbsOptions['gridCenterMethod'] = "coordinate"
        apbsOptions['gridCenter'] = myElec.gcent


    if myElec.lpbe == 1:
        apbsOptions['solveType'] = 'linearized'
    elif myElec.npbe == 1:
        apbsOptions['solveType'] = 'nonlinearized'

    if len(myElec.ion) == 0:
        apbsOptions['mobileIonSpecies[0]'] = None
    else:
        apbsOptions['mobileIonSpecies[1]'] = myElec.ion

    if len(myElec.write) <= 1:
        apbsOptions['format'] = 'dx'
    else:
        apbsOptions['format'] = myElec.write[1]

    apbsOptions['calculationType'] = myElec.method
    apbsOptions['dime'] = myElec.dime
    apbsOptions['pdime'] = myElec.pdime
    apbsOptions['async'] = myElec.async
    apbsOptions['asyncflag'] = myElec.asyncflag
    apbsOptions['nlev'] = myElec.nlev
    apbsOptions['glen'] = myElec.glen
    apbsOptions['coarseGridLength'] = myElec.cglen
    apbsOptions['fineGridLength'] = myElec.fglen
    apbsOptions['molecule'] = myElec.mol
    apbsOptions['boundaryConditions'] = myElec.bcfl
    apbsOptions['biomolecularDielectricConstant'] = myElec.pdie
    apbsOptions['dielectricSolventConstant'] = myElec.sdie
    apbsOptions['biomolecularPointChargeMapMethod'] = myElec.chgm
    apbsOptions['surfaceConstructionResolution'] = myElec.sdens
    apbsOptions['dielectricIonAccessibilityModel'] = myElec.srfm
    apbsOptions['solventRadius'] = myElec.srad
    apbsOptions['surfaceDefSupportSize'] = myElec.swin
    apbsOptions['temperature'] = myElec.temp
    apbsOptions['calculationEnergy'] = myElec.calcenergy
    apbsOptions['calculationForce'] = myElec.calcforce
    apbsOptions['processorMeshOverlap'] = myElec.ofrac
    apbsOptions['writeBiomolecularChargeDistribution'] = False
    apbsOptions['writeElectrostaticPotential'] = True
    apbsOptions['writeMolecularSurfaceSolventAccessibility'] = False
    apbsOptions['writeSplineBasedSolventAccessibility'] = False
    apbsOptions['writeVanDerWaalsSolventAccessibility'] = False
    apbsOptions['writeInflatedVanDerWaalsIonAccessibility'] = False
    apbsOptions['writePotentialLaplacian'] = False
    apbsOptions['writeEnergyDensity'] = False
    apbsOptions['writeMobileIonNumberDensity'] = False
    apbsOptions['writeMobileChargeDensity'] = False
    apbsOptions['writeDielectricMapShift'] = [False,False,False]
    apbsOptions['writeIonAccessibilityKappaMap'] = False

    return apbsOptions

def fieldStorageToDict(form):
    """ Converts the CGI input from the web interface to a dictionary """
    apbsOptions = {'writeCheck':0}

    if form.has_key("writecharge") and form["writecharge"].value != "":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeCharge'] = True
    else:
        apbsOptions['writeCharge'] = False
    
    if form.has_key("writepot") and form["writepot"].value != "":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writePot'] = True
    else:
        apbsOptions['writePot'] = False

    if form.has_key("writesmol") and form["writesmol"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeSmol'] = True
    else:
        apbsOptions['writeSmol'] = False

    if form.has_key("writesspl") and form["writesspl"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeSspl'] = True
    else:
        apbsOptions['writeSspl'] = False

    if form.has_key("writevdw") and form["writevdw"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeVdw'] = True
    else:
        apbsOptions['writeVdw'] = False

    if form.has_key("writeivdw") and form["writeivdw"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeIvdw'] = True
    else:
        apbsOptions['writeIvdw'] = False

    if form.has_key("writelap") and form["writelap"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeLap'] = True
    else:
        apbsOptions['writeLap'] = False

    if form.has_key("writeedens") and form["writeedens"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeEdens'] = True
    else:
        apbsOptions['writeEdens'] = False

    if form.has_key("writendens") and form["writendens"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeNdens'] = True
    else:
        apbsOptions['writeNdens'] = False

    if form.has_key("writeqdens") and form["writeqdens"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeQdens'] = True
    else:
        apbsOptions['writeQdens'] = False

    if form.has_key("writedielx") and form["writedielx"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeDielx'] = True
    else:
        apbsOptions['writeDielx'] = False

    if form.has_key("writediely") and form["writediely"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeDiely'] = True
    else:
        apbsOptions['writeDiely'] = False

    if form.has_key("writedielz") and form["writedielz"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeDielz'] = True
    else:
        apbsOptions['writeDielz'] = False

    if form.has_key("writekappa") and form["writekappa"].value == "on":
        apbsOptions['writeCheck'] += 1
        apbsOptions['writeKappa'] = True
    else:
        apbsOptions['writeKappa'] = False
    
    if apbsOptions['writeCheck'] > 4:
        print "Please select a maximum of four write statements."
        os._exit(99)

    # READ section variables
    apbsOptions['readType'] = "mol"
    apbsOptions['readFormat'] = "pqr"
    apbsOptions['pqrPath'] = ""
    apbsOptions['pqrFileName'] = form['pdb2pqrid'].value+'.pqr'

    #ELEC section variables
    apbsOptions['calcType'] = form["type"].value  

    apbsOptions['dimeNX'] = locale.atoi(form["dimenx"].value)
    apbsOptions['dimeNY'] = locale.atoi(form["dimeny"].value)
    apbsOptions['dimeNZ'] = locale.atoi(form["dimenz"].value)

    apbsOptions['cglenX'] = locale.atof(form["cglenx"].value)
    apbsOptions['cglenY'] = locale.atof(form["cgleny"].value)
    apbsOptions['cglenZ'] = locale.atof(form["cglenz"].value)

    apbsOptions['fglenX'] = locale.atof(form["fglenx"].value)
    apbsOptions['fglenY'] = locale.atof(form["fgleny"].value)
    apbsOptions['fglenZ'] = locale.atof(form["fglenz"].value)

    apbsOptions['glenX'] = locale.atof(form["glenx"].value)
    apbsOptions['glenY'] = locale.atof(form["gleny"].value)
    apbsOptions['glenZ'] = locale.atof(form["glenz"].value)

    if form["cgcent"].value == "mol":
        apbsOptions['coarseGridCenterMethod'] = "molecule"
        apbsOptions['coarseGridCenterMoleculeID'] = locale.atoi(form["cgcentid"].value)

    elif form["cgcent"].value == "coord":
        apbsOptions['coarseGridCenterMethod'] = "coordinate"
        apbsOptions['cgxCent'] = locale.atoi(form["cgxcent"].value)
        apbsOptions['cgyCent'] = locale.atoi(form["cgycent"].value)
        apbsOptions['cgzCent'] = locale.atoi(form["cgzcent"].value)

    if form["fgcent"].value == "mol":
        apbsOptions['fineGridCenterMethod'] = "molecule"
        apbsOptions['fineGridCenterMoleculeID'] = locale.atoi(form["fgcentid"].value)
    elif form["fgcent"].value == "coord":
        apbsOptions['fineGridCenterMethod'] = "coordinate"
        apbsOptions['fgxCent'] = locale.atoi(form["fgxcent"].value)
        apbsOptions['fgyCent'] = locale.atoi(form["fgycent"].value)
        apbsOptions['fgzCent'] = locale.atoi(form["fgzcent"].value)

    if form["gcent"].value == "mol":
        apbsOptions['gridCenterMethod'] = "molecule"
        apbsOptions['gridCenterMoleculeID'] = locale.atoi(form["gcentid"].value)
    elif form["gcent"].value == "coord":
        apbsOptions['gridCenterMethod'] = "coordinate"
        apbsOptions['gxCent'] = locale.atoi(form["gxcent"].value)
        apbsOptions['gyCent'] = locale.atoi(form["gycent"].value)
        apbsOptions['gzCent'] = locale.atoi(form["gzcent"].value)


    apbsOptions['mol'] = locale.atoi(form["mol"].value)
    apbsOptions['solveType'] = form["solvetype"].value
    apbsOptions['boundaryConditions'] = form["bcfl"].value
    apbsOptions['biomolecularDielectricConstant'] = locale.atof(form["pdie"].value)
    apbsOptions['dielectricSolventConstant'] = locale.atof(form["sdie"].value)
    apbsOptions['dielectricIonAccessibilityModel'] = form["srfm"].value
    apbsOptions['biomolecularPointChargeMapMethod'] = form["chgm"].value
    apbsOptions['surfaceConstructionResolution'] = locale.atof(form["sdens"].value)
    apbsOptions['solventRadius'] = locale.atof(form["srad"].value)    
    apbsOptions['surfaceDefSupportSize'] = locale.atof(form["swin"].value)
    apbsOptions['temperature'] = locale.atof(form["temp"].value)
    apbsOptions['calcEnergy'] = form["calcenergy"].value
    apbsOptions['calcForce'] = form["calcforce"].value

    for i in range(0,3):
        if form['charge%i' % i].value != "":
            apbsOptions['mobileIonSpeciesCharge'] = locale.atoi(form['charge%i' % i].value)
        if form['conc%i' % i].value != "":
            apbsOptions['mobileIonSpeciesConcentration'] = locale.atof(form['conc%i' % i].value)
        if form['radius%i' % i].value != "":
            apbsOptions['mobileIonSpeciesRadius'] = locale.atof(form['radius%i' % i].value)
    apbsOptions['writeFormat'] = form["writeformat"].value
    #apbsOptions['writeStem'] = apbsOptions['pqrFileName'][:-4]
    apbsOptions['writeStem'] = form["pdb2pqrid"].value


    return apbsOptions


def pqrFileCreator(apbsOptions):
    """
        Creates a pqr file, using the data from the form
    """
    apbsOptions['tmpDirName'] = "%s%s%s/" % (INSTALLDIR, TMPDIR, apbsOptions['writeStem'])
    try:
        os.makedirs(apbsOptions['tmpDirName'])
    except OSError, err:
        if err.errno == errno.EEXIST:
            if os.path.isdir(apbsOptions['tmpDirName']):
                # print "Error (tmp directory already exists) - please try again"
                pass
            else:
                print "Error (file exists where tmp dir should be) - please try again"
                raise
        else:
            raise

    apbsOptions['tempFile'] = "apbsinput.in"
    apbsOptions['tab'] = "    " # 4 spaces - used for writing to file
    input = open('%s/tmp/%s/%s' % (INSTALLDIR, apbsOptions['writeStem'], apbsOptions['tempFile']), 'w')
    

    # writing READ section to file
    input.write('read\n')
    input.write('%s%s %s %s%s' % (apbsOptions['tab'], apbsOptions['readType'], apbsOptions['readFormat'], apbsOptions['pqrPath'], apbsOptions['pqrFileName']))
    input.write('\nend\n')

    # writing ELEC section to file
    input.write('elec\n')
    input.write('%s%s\n' % (apbsOptions['tab'], apbsOptions['calcType']))
    if apbsOptions['calcType']!="fe-manual":
        input.write('%sdime %d %d %d\n' % (apbsOptions['tab'], apbsOptions['dimeNX'], apbsOptions['dimeNY'], apbsOptions['dimeNZ']))
    if apbsOptions['calcType'] == "mg-para":
        input.write('%sdime %d %d %d\n' % (apbsOptions['tab'], apbsOptions['pdime'][0], apbsOptions['pdime'][1], apbsOptions['pdime'][2]))

    if apbsOptions['calcType'] == "mg-manual":
        input.write('%sglen %g %g %g\n' % (apbsOptions['tab'], apbsOptions['glenX'], apbsOptions['glenY'], apbsOptions['glenZ']))
    if apbsOptions['calcType'] in ['mg-auto','mg-para','mg-dummy']:
        input.write('%scglen %g %g %g\n' % (apbsOptions['tab'], apbsOptions['cglenX'], apbsOptions['cglenY'], apbsOptions['cglenZ']))
    if apbsOptions['calcType'] in ['mg-auto','mg-para']:
        input.write('%sfglen %g %g %g\n' % (apbsOptions['tab'], apbsOptions['fglenX'], apbsOptions['fglenY'], apbsOptions['fglenZ']))

        if apbsOptions['coarseGridCenterMethod']=='molecule':
            input.write('%scgcent mol %d\n' % (apbsOptions['tab'], apbsOptions['coarseGridCenterMoleculeID'] ))
        elif apbsOptions['coarseGridCenterMethod']=='coordinate':
            input.write('%scgcent %d %d %d\n' % (apbsOptions['tab'], apbsOptions['cgxCent'], apbsOptions['cgyCent'], apbsOptions['cgzCent']))

        if apbsOptions['fineGridCenterMethod']=='molecule':
            input.write('%sfgcent mol %d\n' % (apbsOptions['tab'], apbsOptions['fineGridCenterMoleculeID']))
        elif apbsOptions['fineGridCenterMethod']=='coordinate':
            input.write('%sfgcent %d %d %d\n' % (apbsOptions['tab'], apbsOptions['fgxCent'], apbsOptions['fgyCent'], apbsOptions['fgzCent']))

    if apbsOptions['calcType'] in ['mg-manual','mg-dummy']:
        if apbsOptions['gridCenterMethod']=='molecule':
            input.write('%sgcent mol %d\n' % (apbsOptions['tab'], apbsOptions['gridCenterMoleculeID'] ))
        elif apbsOptions['gridCenterMethod']=='coordinate':
            input.write('%sgcent %d %d %d\n' % (apbsOptions['tab'], apbsOptions['gxCent'], apbsOptions['gyCent'], apbsOptions['gzCent']))

    input.write('%smol %d\n' % (apbsOptions['tab'], apbsOptions['mol']))
    input.write('%s%s\n' % (apbsOptions['tab'], apbsOptions['solveType']))
    input.write('%sbcfl %s\n' % (apbsOptions['tab'], apbsOptions['boundaryConditions']))
    input.write('%spdie %g\n' % (apbsOptions['tab'], apbsOptions['biomolecularDielectricConstant']))
    input.write('%ssdie %g\n' % (apbsOptions['tab'], apbsOptions['dielectricSolventConstant']))
    input.write('%ssrfm %s\n' % (apbsOptions['tab'], apbsOptions['dielectricIonAccessibilityModel']))
    input.write('%schgm %s\n' % (apbsOptions['tab'], apbsOptions['biomolecularPointChargeMapMethod']))
    input.write('%ssdens %g\n' % (apbsOptions['tab'], apbsOptions['surfaceConstructionResolution']))
    input.write('%ssrad %g\n' % (apbsOptions['tab'], apbsOptions['solventRadius']))
    input.write('%sswin %g\n' % (apbsOptions['tab'], apbsOptions['surfaceDefSupportSize']))
    input.write('%stemp %g\n' % (apbsOptions['tab'], apbsOptions['temperature']))
    input.write('%scalcenergy %s\n' % (apbsOptions['tab'], apbsOptions['calcEnergy']))
    input.write('%scalcforce %s\n' % (apbsOptions['tab'], apbsOptions['calcForce']))
    if apbsOptions.has_key('mobileIonSpeciesCharge') and apbsOptions.has_key('conc') and apbsOptions.has_key('radius'):
        input.write('%sion %d %g %g\n' % (apbsOptions['tab'], apbsOptions['mobileIonSpeciesCharge'], apbsOptions['mobileIonSpeciesConcentration'], apbsOptions['mobileIonSpeciesRadius']))


    if apbsOptions['writeCharge']:
        input.write('%swrite charge %s %s-charge\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))
    
    if apbsOptions['writePot']:
        input.write('%swrite pot %s %s-pot\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeSmol']:
        input.write('%swrite smol %s %s-smol\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeSspl']:
        input.write('%swrite sspl %s %s-sspl\n' % (apbsOptions['tab'], apbsOptions['writeFormat'],  apbsOptions['writeStem']))

    if apbsOptions['writeVdw']:
        input.write('%swrite vdw %s %s-vdw\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeIvdw']:
        input.write('%swrite ivdw %s %s-ivdw\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeLap']:
        input.write('%swrite lap %s %s-lap\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeEdens']:
        input.write('%swrite edens %s %s-edens\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeNdens']:
        input.write('%swrite ndens %s %s-ndens\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeQdens']:
        input.write('%swrite qdens %s %s-qdens\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeDielx']:
        input.write('%swrite dielx %s %s-dielx\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeDiely']:
        input.write('%swrite diely %s %s-diely\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeDielz']:
        input.write('%swrite dielz %s %s-dielz\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    if apbsOptions['writeKappa']:
        input.write('%swrite kappa %s %s-kappa\n' % (apbsOptions['tab'], apbsOptions['writeFormat'], apbsOptions['writeStem']))

    input.write('end\n')
    input.write('quit')
    input.close()





def convertOpalToLocal(jobid,pdb2pqrOpalJobID):
    """
        takes a remote Opal page and saves the files to a local directory
    """
    appLocator = AppServiceLocator()
    resp = appLocator.getAppServicePort(PDB2PQR_OPAL_URL).getOutputs(getOutputsRequest(pdb2pqrOpalJobID))

    # returns the variable to prepend the files
    #sys.path.append('/home/samir/public_html/pdb2pqr/src') # HARDCODED
    #from src import server
    #logTime = setID(time.time())
    #os.makedirs('%s%s%s' % (INSTALLDIR, TMPDIR, logTime))
    for file in resp._outputFile:
        fileName = file._name
        if fileName!="Standard Output" and fileName!="Standard Error":
            if fileName.rfind('-') != -1:
                fileName = jobid+fileName[fileName.rfind('-'):]
            elif fileName.rfind('.') != -1:
                fileName = jobid+fileName[fileName.rfind('.'):]
            urllib.urlretrieve(file._url, '%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, fileName)) # HARDCODED

def redirector(logTime):
    if (str(logTime) != "False") and (str(logTime) != "notenoughmem"):
        starttimefile = open('%s%s%s/apbs_start_time' % (INSTALLDIR, TMPDIR, logTime), 'w')
        starttimefile.write(str(time.time()))
        starttimefile.close()

    string = ""
    string+='<html> <head>'
    # status is passed to querystatus.cgi
    string+='<meta http-equiv=\"refresh\" content=\"0;url=querystatus.cgi?jobid=%s&calctype=apbs\"/></head></html>' % str(logTime)
    return string

def mainInput() :
    """
        Main function
    """
    global have_opal
    file = stdout
    file.write("Content-type: text/html; charset=utf-8\n\n") 
    cgitb.enable()

    # Check cgi.FieldStorage() for checkbox indicating whether we were invoked frmo the form or from the URL
    form = cgi.FieldStorage()
    
    firstRun = True
    pdb2pqrOpalChecked = False
    pdb2pqrChecked = False

    if form.has_key("jobid"): # means it's not the first run
        firstRun = False
        if HAVE_PDB2PQR_OPAL=="1": 
            pdb2pqrOpalJobIDFile = open('%s%s%s/pdb2pqr_opal_job_id' % (INSTALLDIR, TMPDIR, form["jobid"].value))
            pdb2pqrOpalJobID = pdb2pqrOpalJobIDFile.read()
            pdb2pqrOpalJobIDFile.close()
            pdb2pqrID = form["jobid"].value
            logTime = form["jobid"].value
            convertOpalToLocal(form["jobid"].value, pdb2pqrOpalJobID)
            pdb2pqrOpalChecked = True
            pdb2pqrChecked = False
        else:
            pdb2pqrID = form["jobid"].value
            logTime = pdb2pqrID
            pdb2pqrChecked = True
            pdb2pqrOpalChecked = False

    if form.has_key("hiddencheck"): # means this time apbs must be run
        firstRun = False
        if form["hiddencheck"].value == "local":
            typeOfRun="local"
        elif form["hiddencheck"].value == "opal":
            typeOfRun="opal"


    #redirects to pdb2pqr input page
    if firstRun:
        pdb2pqrLocation = '../pdb2pqr/html/server.html'
        print '<html>'
        print '<head>'
        print '<meta http-equiv=\"refresh\" content=\"0;url=%s\"/>' % pdb2pqrLocation
        print '</head>'
        print '<body>'
        print '</body>'
        print '</html>'


    #generates web interface and displays it
    elif pdb2pqrChecked:
        initVars = unpickleVars(pdb2pqrID)

        generateForm(file, initVars, pdb2pqrID, "local")

    elif pdb2pqrOpalChecked:
        initVars = unpickleVars(pdb2pqrID)

        generateForm(file, initVars, pdb2pqrID, "opal")

    #runs apbs
    else:
        # logTime stores the prefix of the names for all the data files for a run
        logTime = form["pdb2pqrid"].value

        tempPage = "results.html"   

        apbsOptions = fieldStorageToDict(form)
        pqrFileCreator(apbsOptions)
        if APBS_OPAL_URL == "":
            have_opal = False
        else:
            have_opal = True

        aoFile = open('%s%s%s/%s-ao' % (INSTALLDIR, TMPDIR, logTime, logTime),'w')
        pickle.dump(apbsOptions, aoFile)
        aoFile.close()

        if have_opal:
            apbsOpalJobID = apbsOpalExec(logTime, form, apbsOptions)

            # if the version number doesn't match, apbsOpalExec returns False
            if(str(apbsOpalJobID) == 'False'):
                print redirector(False)

            # Check if not enough memory
            elif(str(apbsOpalJobID) == 'notenoughmem'):
                print redirector('notenoughmem')
            else:
                print redirector(logTime)

            apbsOpalJobIDFile = open('%s%s%s/apbs_opal_job_id' % (INSTALLDIR, TMPDIR, logTime),'w')
            apbsOpalJobIDFile.write(apbsOpalJobID)
            apbsOpalJobIDFile.close()
        else:
            apbsExec(logTime, form, apbsOptions)


if __name__ == "__main__" and os.environ.has_key("REQUEST_METHOD"):
    """ Determine if called from command line or CGI """

    if APBS_OPAL_URL!="" or HAVE_PDB2PQR_OPAL=="1":
        have_opal = True
        from AppService_client import queryStatusRequest
        from AppService_client import AppServiceLocator, queryStatusRequest, getOutputsRequest
    else:
        have_opal = False
    mainInput()
