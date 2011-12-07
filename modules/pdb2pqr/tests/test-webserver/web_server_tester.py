"""
    Web server testing code for PDB2PQR

    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

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

    ----------------------------
"""

__date__  = "26 January 2010"
__author__ = "Samir Unni, Yong Huang"

import urllib, urllib2, time, StringIO
import os
import random
from sys import argv

rand = int(random.random()*3) + 1     # Creating a random integer in {1,2,3}

def launchJob(pdbID):
    url = 'http://pdb2pqr-%s/pdb2pqr/pdb2pqr.cgi' % (rand)
    values = {'TEXTCONTROL' : 'ID',
              'PDBID' : pdbID,
              'FF' : 'amber',
              'FFOUT' : 'internal',
              'DEBUMP' : 'CHECKED',
              'OPT' : 'CHECKED'}
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)

    for i in range(0,2):
        response.readline()
    url = response.readline()[45:-3]
    return url

def checkJob(url):
    response = open(urllib.urlretrieve(url)[0])
    while True:
        status = response.readline()
        if status[0:7]=="Message":
            status = status[9:-7].strip()
            break
        elif status=='':
            break

    return status

def main():

    #---User configurable values---
    threads = 100 # number of simultaneous calculations
    timeout = 1000 # timeout in seconds (before a job is considered failed)
    sleeptime = 30 # how long to wait between checks in seconds
    #------------------------

    currentdir = os.getcwd()
    webtestdir = currentdir + '/tests/test-webserver/'
    listFile = open(webtestdir+'lst','r')
    pdbIDs = listFile.read().strip().split('\n')
    listFile.close()
    currentID = 0
    currentJobs = []
    completedJobs = []
    failedJobs = []
    timeout = timeout/sleeptime

    while(currentID<len(pdbIDs) or len(currentJobs)>0):
        while(len(currentJobs)<threads and currentID<len(pdbIDs)):
            pdbID = pdbIDs[currentID].strip()
            print "launching %s" % pdbID
            currentJobs.append([pdbID,launchJob(pdbID),0])
            print "%i threads" % threads
            currentID = currentID+1
        print "sleeping"
        time.sleep(sleeptime)
        jobsToDelete = []
        for i in range(0,len(currentJobs)):
            print currentJobs[i]
            print "checking %s" % currentJobs[i][0]
            result = checkJob(currentJobs[i][1])
            print "%s:  %s" % (currentJobs[i][0], result)
            if(result=="complete" or result=="failed"):
                completedJobs.append([currentJobs[i][0], result])
                jobsToDelete.append(i)
            elif currentJobs[i][2]>timeout:
                print "%s: failed" % currentJobs[i][0]
                completedJobs.append([currentJobs[i][0], "failed"])
                jobsToDelete.append(i)
            else:
                currentJobs[i][2] = currentJobs[i][2]+1
        print "jobsToDelete: %s" % jobsToDelete
        jobsToDelete.reverse()
        for i in jobsToDelete:
            print "i = %i" % i
            print "deleting %s" % currentJobs[i][0]
            del currentJobs[i]

    for i in range(0,len(completedJobs)):
        if(completedJobs[i][1]=="failed"):
            failedJobs.append(completedJobs[i][0])

    if len(failedJobs) == 0:
        print "\nWeb test finished with no failed jobs.\n"
    else:
        print "\nThese jobs failed: %s.\n" % (failedJobs) 

if __name__ == "__main__":
    main()
