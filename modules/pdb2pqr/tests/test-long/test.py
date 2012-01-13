"""
    PDB2PQR test on a long list of PDB IDs

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

__date__   = "14 May 2009"
__author__ = "Yong Huang"

import os
import sys
import getopt
import string

def runtest(argv):
    
    count = 0
    currentdir = os.getcwd()
    longtestdir = currentdir + '/tests/test-long/'
    f=open(longtestdir+'lst')
    line=f.readlines()
    defaultnum = len(line)   # Default numbers of PDB tests in test-long

    options={"testnum": defaultnum}
 
    try: opts, args = getopt.getopt(sys.argv[1:], 'n', ['testnum='])
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)

    for o,a in opts:
        if o in ("-n", "--testnum"):
            options["testnum"] = int(a)
        if options["testnum"] >= defaultnum or options["testnum"] <= 0:
            raise ValueError, "TESTNUM must be an integer between 1 and %s!\n" % (defaultnum - 1)

    for element in line:
        
        inname = element.strip()
        outname = os.path.join(longtestdir+'out', inname+'.pqr')
        errname = os.path.join(longtestdir+'out', inname+'.err')
        
        os.system("python pdb2pqr.py --ff=parse %s %s 1>/dev/null 2>%s" % (inname, outname, errname))
        
        count += 1
        
        print "Finished testing for PDB: %s, number of PDBs tested: %s." % (inname, count)

        if 0 < options["testnum"] < defaultnum and count == options["testnum"]:
            break

    print "\nLong test finished, please check tests/test-long/out/ directory for output files.\n"

if __name__ == "__main__":
    import os
    pwd =os.getcwd()

    runtest(sys.argv) 
