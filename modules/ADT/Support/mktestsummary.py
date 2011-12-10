import os, sys
import time
import stat
import string

mglpacks=[ "DejaVu", "MolKit", "ViewerFramework",  "Pmv", "AutoDockTools", "PyBabel",
           "symserv", "mglutil", "NetworkEditor", "Volume", "Vision","PyAutoDock", "Support"]                              
    
def getTestsDirs(file_list, directory, names):
    if os.path.basename(directory) == "Tests":
        file_list.append(os.path.dirname(directory))
        
def getTestSummaryResults(testroot, testreportdir):
    """Parse the test report files created by tester.
    testroot - directory containing the MGLTools packages ;
    testreportdir - directory containing testreports """
    #print "testreportdir:", testreportdir
    assert os.path.exists(testreportdir)
    assert os.path.exists(testroot)   
    # list of directories containing Tests ( created when the tests ran last for this platform)

    cwd = os.getcwd()
    
    testdirs = []
    os.chdir(testroot)
    for p in mglpacks :
        os.path.walk(p, getTestsDirs, testdirs)
    #print  "testdirs" ,testdirs
    os.chdir(cwd)
    reports = []
    for p in testdirs:
        reports.append(p.replace(os.sep,'.'))

    # a dictionary containing the names of test reports
    rantests = {}
    # a dictionary containing the names of test reports (with -s option)
    rantestsS = {}
    # a list of test modules for which we did not get a report
    noreports= []

    modules = []
    # string containing time when the directory was last modified
    reporttime = time.asctime(time.localtime(os.stat(testreportdir)[stat.ST_MTIME])) 

    # remove time (##:##:##)   from reporttime string (leave only the date)
    rd = string.split(reporttime)
    rd.pop(-2)
    reportdate = "%s %s %s %s" % tuple(rd)


    import sys
    for r in reports:
        f1 = os.path.join(testroot, testreportdir,
                          "report-" + r + sys.platform + ".txt")
        f2 = os.path.join(testroot, testreportdir,
                          "sReport-"+ r + sys.platform + ".txt")

        if os.path.exists(f1):
            rantests[f1] = r 
        elif os.path.exists(f2):
            rantestsS[f2] = r
        else:
            noreports.append(r)

    numtests = 0
    numerrors = 0
    results = {}
    for file in rantestsS.keys():
        testerrors = 0
        f= open(file, "r")
        st =f.readline()
        l = string.split(st)
        if len(l)== 5:
            if l[0] == "RAN":
                test = rantestsS[file]
                ntests = float(l[1])
                numtests = numtests + ntests
                modules.append(test)
                results[test]={}
                results[test]['ntests']= ntests
        txt = f.readlines()
        f.close()
        testname = None
        errlist = []
        errtxt = ""
        for line in txt:
            if string.find(line, "Tests.test_")> 0:
                if string.find(line, "ERROR:")== 0 or string.find(line, "FAIL:")==0:
                    errlist.append(line)
                    continue
                if len(errlist) >0:
                    errtxt = errtxt + "\n" + testname
                    for err in errlist:
                        errtxt = errtxt + err
                    testerrors = testerrors +len(errlist)
                    errlist = []
                testname =  line
                continue
            elif string.find(line, "ERROR:")== 0 or string.find(line, "FAIL:")==0:
                errlist.append(line)
        #this is for the case when error is in the last test:       
        if len(errlist) >0:
            errtxt = errtxt + "\n" + testname
            for err in errlist:
                errtxt = errtxt + err
            testerrors = testerrors +len(errlist)
            errlist = []

        if testerrors:
            numerrors = numerrors + testerrors
            results[test]['failed'] = testerrors
            results[test]['errtxt'] = errtxt

    for file in rantests.keys():
        testerrors = 0
        f= open(file, "r")
        st =f.readline()
        l = string.split(st)
        if len(l)== 5:
            if l[0] == "RAN":
                test = rantests[file]
                modules.append(test)
                ntests = float(l[1])
                numtests = numtests + ntests
                results[test]={}
                results[test]['ntests']= ntests

        txt = f.readlines()
        f.close()
        testname = None
        errlist = []
        errtxt = ""
        for line in txt:
            if string.find(line, "FAILED (")==0:
                if len(errlist) >0:
                    for err in errlist:
                        errtxt = errtxt + err
                    testerrors = testerrors +len(errlist)
                    errlist = []
                continue
            elif string.find(line, "ERROR:")== 0 or string.find(line, "FAIL:")==0:
                errlist.append(line)
        if testerrors:
            numerrors = numerrors + testerrors
            results[test]['failed'] = testerrors
            results[test]['errtxt'] = errtxt

    perc = numerrors/numtests*100.0
    if len(noreports):
        endstr = "No test report for the following module(s) (tests possibly crashed): \n"
        for f in noreports:
            endstr = endstr+ "%s \n" % f
        results['tail']= endstr   
    results['reportdate'] = reportdate
    results['numtests'] = numtests
    results['numerrors'] = numerrors
    results['perc'] = perc

    return results, modules



def mkSummaryHTML(testroot, testreportdir, filename = None):
    """ create a simple html file with the tests summary"""
    
    results, modules = getTestSummaryResults(testroot, testreportdir)
    if not filename:
        if sys.platform == "darwin":
            if os.uname()[-1] == 'Power Macintosh':
                prt = "powerpc"
            else:
                prt = "i386"
            filename = "testsummary"+sys.platform+ prt+".html"
        else:
            filename = "testsummary"+sys.platform+".html"
        
    resfile = open(filename, "w")
    
    txt = """
    <!-- <html>
    <head>
    <title>TestSummary</title>
    </head>
    <body bgcolor=#ffffff text=#000000> -->
    <html>
    <BODY text=#000000 bgcolor=#FFFFFF >
    <font color=#7c0027 size=+1>Summary of Tests Results ran on %s.</font>
    <br> """ % results['reportdate']
    resfile.write(txt)
    for mod in modules:
        ntests = results[mod]['ntests']
        if ntests == 0:
            continue
        #resfile.write("<b>%s</b><br> \n" % mod)
        resfile.write("<b>%s &nbsp; &nbsp; </b>\n" % mod+"   ")
        
        if results[mod].has_key("failed"):
            nfailed = results[mod]["failed"]
            #txt = "<pre>\nRan <b>%d</b> tests: number of failed tests <b>%d</b> \nThe following tests failed: \n" % (ntests, nfailed)
            txt = "Ran <b>%d</b> tests; <b>%d</b> test(s) failed:\n" % (ntests, nfailed)
            resfile.write(txt)
            resfile.write("<pre>")
            resfile.write(results[mod]['errtxt']+ "</pre> \n" )
            
            #resfile.write("</pre>\n-----------------------------------------------------------<br>\n")
        else:
            resfile.write("Ran <b>%d</b> tests: OK <br>\n" % ntests)
            #resfile.write("-----------------------------------------------------------<br>\n")

    resfile.write("===========================================================<br>\n")
    resfile.write("<b>Overall Results(%s): ran %d tests, failed %d (%.2f %s) </b><br> \n" %
                  (results['reportdate'], results['numtests'], results['numerrors'], results['perc'], "%") )
    if results.has_key('tail'):
        resfile.write("<pre>\n%s</pre>" % results['tail'])
    
    resfile.close()
    

if __name__ == "__main__":

    args = sys.argv[1:]
    #print "args" ,args
    if len(args) == 2:
        mkSummaryHTML(args[0], args[1])
    elif len(args) == 3:
        mkSummaryHTML(args[0], args[1], args[2])
    else:
        print "mkreport requires 2 or 3 arguments"

