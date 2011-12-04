""" psize class

    Get dimensions and other information from a PQR file.

    Originally written by Dave Sept
    Additional APBS-specific features added by Nathan Baker
    Ported to Python/Psize class by Todd Dolinsky and subsequently
    hacked by Nathan Baker

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

# User - Definable Variables: Default values

# CFAC = 1.7                  # Factor by which to expand mol dims to
                              # get coarse grid dims
# FADD = 20                   # Amount to add to mol dims to get fine
                              # grid dims
# SPACE = 0.50                # Desired fine mesh resolution
# GMEMFAC = 200               # Number of bytes per grid point required 
                              # for sequential MG calculation 
# GMEMCEIL = 400              # Max MB allowed for sequential MG
                              # calculation.  Adjust this to force the
                              # script to perform faster calculations (which
                              # require more parallelism).
# OFAC = 0.1                  # Overlap factor between mesh partitions
# REDFAC = 0.25               # The maximum factor by which a domain
                              # dimension can be reduced during focusing
# TFAC_ALPHA = 9e-5           # Number of sec/unkown for setup/solve on 667
                              # MHz EV67 Alpha CPU -- VERY ROUGH ESTIMATE
# TFAC_XEON  = 3e-4           # Number of sec/unkown for setup/solve on 500
                              # MHz PIII Xeon CPU -- VERY ROUGH ESTIMATE
# TFAC_SPARC  = 5e-4          # Number of sec/unkown for setup/solve on 400
                              # MHz UltraSPARC II CPU -- VERY ROUGH ESTIMATE

import string, sys
from sys import stdout
from math import log

class Psize:
    def __init__(self):
        self.constants = {"CFAC":1.7, "FADD":20, "SPACE":0.50, "GMEMFAC":200, "GMEMCEIL":400, "OFAC":0.1, "REDFAC":0.25, "TFAC_ALPHA":9e-5, "TFAC_XEON":3e-4, "TFAC_SPARC": 5e-4}
        self.minlen = [360.0, 360.0, 360.0]
        self.maxlen = [0.0, 0.0, 0.0]
        self.q = 0.0
        self.gotatom = 0
        self.gothet = 0
        self.olen = [0.0, 0.0, 0.0]
        self.cen = [0.0, 0.0, 0.0]
        self.clen = [0.0, 0.0, 0.0]
        self.flen = [0.0, 0.0, 0.0]
        self.n = [0, 0, 0]
        self.np = [0.0, 0.0, 0.0]
        self.nsmall = [0,0,0]
        self.nfocus = 0

    def parseString(self, structure):
        """ Parse the input structure as a string in PDB or PQR format """
        lines = string.split(structure, "\n")
        self.parseLines(lines)

    def parseInput(self, filename):
        """ Parse input structure file in PDB or PQR format """
        file = open(filename, "r")
        self.parseLines(file.readlines())

    def parseLines(self, lines):
        """ Parse the lines """
        for line in lines:
            if string.find(line,"ATOM") == 0:
                subline = string.replace(line[30:], "-", " -")
                words = string.split(subline)
                if len(words) < 4:    
                    #sys.stderr.write("Can't parse following line:\n")
                    #sys.stderr.write("%s\n" % line)
                    #sys.exit(2)
                    continue
                self.gotatom = self.gotatom + 1
                self.q = self.q + float(words[3])
                rad = float(words[4])
                center = []
                for word in words[0:3]:
                    center.append(float(word))
                for i in range(3):
                    self.minlen[i] = min(center[i]-rad, self.minlen[i])
                    self.maxlen[i] = max(center[i]+rad, self.maxlen[i])
            elif string.find(line, "HETATM") == 0:
                self.gothet = self.gothet + 1
    
    def setConstant(self, name, value):
        """ Set a constant to a value; returns 0 if constant not found """
        try:
            self.constants[name] = value
            return 1
        except KeyError:
            return 0

    def getConstant(self, name):
        """ Get a constant value; raises KeyError if constant not found """
        return self.constants[name]

    def setLength(self, maxlen, minlen):
        """ Compute molecule dimensions """
        for i in range(3):
            self.olen[i] = maxlen[i] - minlen[i]
            if self.olen[i] < 0.1:
                self.olen[i] = 0.1
        return self.olen


    def setCoarseGridDims(self, olen):
        """ Compute coarse mesh dimensions """
        for i in range(3):
            self.clen[i] = self.constants["CFAC"] * olen[i]
        return self.clen

    def setFineGridDims(self, olen, clen):
        """ Compute fine mesh dimensions """
        for i in range(3):
            self.flen[i] = olen[i] + self.constants["FADD"]
            if self.flen[i] > clen[i]:
                #str = "WARNING: Fine length (%.2f) cannot be larger than coarse length (%.2f)\n" % (self.flen[i], clen[i])
                #str = str + "         Setting fine grid length equal to coarse grid length\n"
                #stdout.write(str)
                self.flen[i] = clen[i]
        return self.flen

    def setCenter(self, maxlen, minlen):
        """ Compute molecule center """
        for i in range(3):
            self.cen[i] = (maxlen[i] + minlen[i]) / 2
        return self.cen


    def setFineGridPoints(self, flen):
        """ Compute mesh grid points, assuming 4 levels in MG hierarchy """
        tn = [0,0,0]
        for i in range(3):
            tn[i] = int(flen[i]/self.constants["SPACE"] + 0.5)
            self.n[i] = 32*(int((tn[i] - 1) / 32 + 0.5)) + 1
            if self.n[i] < 33:
                self.n[i] = 33
        return self.n

    def setSmallest(self, n):
        """ Compute parallel division in case memory requirement above ceiling
        Find the smallest dimension and see if the number of grid points in
        that dimension will fit below the memory ceiling
        Reduce nsmall until an nsmall^3 domain will fit into memory """
        nsmall = []
        for i in range(3):
            nsmall.append(n[i])
        while 1:
            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            if nsmem < self.constants["GMEMCEIL"]: break
            else:
                i = nsmall.index(max(nsmall))
                nsmall[i] = 32 * ((nsmall[i] - 1)/32 - 1) + 1
                if nsmall <= 0:
                    stdout.write("You picked a memory ceiling that is too small\n")
                    sys.exit(0)        

        self.nsmall = nsmall
        return nsmall

    def setProcGrid(self, n, nsmall):
        """ Calculate the number of processors required to span each 
        dimension """

        zofac = 1 + 2 * self.constants["OFAC"]
        for i in range(3):
            self.np[i] = n[i]/float(nsmall[i])
            if self.np[i] > 1: self.np[i] = int(zofac*n[1]/nsmall[i] + 1.0)
        return self.np
                                                
    def setFocus(self, flen, np, clen):
        """ Calculate the number of levels of focusing required for each
        processor subdomain """

        nfoc = [0,0,0]
        for i in range(3):
            nfoc[i] = int(log((flen[i]/np[i])/clen[i])/log(self.constants["REDFAC"]) + 1.0)
        nfocus = nfoc[0]
        if nfoc[1] > nfocus: nfocus = nfoc[1]
        if nfoc[2] > nfocus: nfocus = nfoc[2]
        if nfocus > 0: nfocus = nfocus + 1
        self.nfocus = nfocus

    def setAll(self):
        """ Set up all of the things calculated individually above """
        maxlen = self.getMax()
        minlen = self.getMin()
        self.setLength(maxlen, minlen)
        olen = self.getLength()
        
        self.setCoarseGridDims(olen)
        clen = self.getCoarseGridDims()        
        
        self.setFineGridDims(olen, clen)
        flen = self.getFineGridDims()
        
        self.setCenter(maxlen, minlen)
        cen = self.getCenter()
        
        self.setFineGridPoints(flen)
        n = self.getFineGridPoints()
        
        self.setSmallest(n)
        nsmall = self.getSmallest()
        
        self.setProcGrid(n, nsmall)
        np = self.getProcGrid()
        
        self.setFocus(flen, np, clen)
        nfocus = self.getFocus()
        
    def getMax(self): return self.maxlen
    def getMin(self): return self.minlen
    def getCharge(self): return self.q
    def getLength(self): return self.olen
    def getCoarseGridDims(self): return self.clen
    def getFineGridDims(self): return self.flen
    def getCenter(self): return self.cen
    def getFineGridPoints(self): return self.n
    def getSmallest(self): return self.nsmall
    def getProcGrid(self): return self.np
    def getFocus(self): return self.nfocus

    def runPsize(self, filename):
        """ Parse input PQR file and set parameters """
        self.parseInput(filename)
        self.setAll()
    
    def printResults(self):
        """ Return a string with the formatted results """

        str = "\n"
        
        if self.gotatom > 0:

            maxlen = self.getMax()
            minlen = self.getMin()
            q = self.getCharge()
            olen = self.getLength()
            clen = self.getCoarseGridDims()        
            flen = self.getFineGridDims()
            cen = self.getCenter()
            n = self.getFineGridPoints()
            nsmall = self.getSmallest()
            np = self.getProcGrid()
            nfocus = self.getFocus()

            # Compute memory requirements

            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            gmem = 200.0 * n[0] * n[1] * n[2] / 1024 / 1024
            
            # Calculate VERY ROUGH wall clock times

            tsolve_alpha = nfocus*nsmall[0]*nsmall[1]*nsmall[2]*self.constants["TFAC_ALPHA"];
            tsolve_xeon = nfocus*nsmall[0]*nsmall[1]*nsmall[2]*self.constants["TFAC_XEON"];
            tsolve_sparc = nfocus*nsmall[0]*nsmall[1]*nsmall[2]*self.constants["TFAC_SPARC"];

            # Print the calculated entries
            str = str + "################# MOLECULE INFO ####################\n"
            str = str + "Number of ATOM entries = %i\n" % self.gotatom
            str = str + "Number of HETATM entries (ignored) = %i\n" % self.gothet
            str = str + "Total charge = %.3f e\n" % q
            str = str + "Dimensions = %.3f x %.3f x %.3f A\n" % (olen[0], olen[1], olen[2])
            str = str + "Center = %.3f x %.3f x %.3f A\n" % (cen[0], cen[1], cen[2])
            str = str + "Lower corner = %.3f x %.3f x %.3f A\n" % (minlen[0], minlen[1], minlen[2])
            str = str + "Upper corner = %.3f x %.3f x %.3f A\n" % (maxlen[0], maxlen[1], maxlen[2])

            str = str + "\n"
            str = str + "############## GENERAL CALCULATION INFO #############\n"
            str = str + "Coarse grid dims = %.3f x %.3f x %.3f A\n" % (clen[0],
clen[1], clen[2])
            str = str + "Fine grid dims = %.3f x %.3f x %.3f A\n" % (flen[0], flen[1], flen[2])
            str = str + "Num. fine grid pts. = %i x %i x %i\n" % (n[0], n[1], n[2])

            if gmem > self.constants["GMEMCEIL"]:
                str = str + "Parallel solve required (%.3f MB > %.3f MB)\n" % (gmem, self.constants["GMEMCEIL"])
                str = str + "Total processors required = %i\n" % (np[0]*np[1]*np[2])
                str = str + "Proc. grid = %i x %i x %i\n" % (np[0], np[1], np[2])
                str = str + "Grid pts. on each proc. = %i x %i x %i\n" % (nsmall[0], nsmall[1], nsmall[2])
                xglob = np[0]*round(nsmall[0]/(1 + 2*self.constants["OFAC"]) - .001)
                yglob = np[1]*round(nsmall[1]/(1 + 2*self.constants["OFAC"]) - .001)
                zglob = np[2]*round(nsmall[2]/(1 + 2*self.constants["OFAC"]) - .001)
                if np[0] == 1: xglob = nsmall[0]
                if np[1] == 1: yglob = nsmall[1]
                if np[2] == 1: zglob = nsmall[2]
                str = str + "Fine mesh spacing = %g x %g x %g A\n" % (flen[0]/(xglob-1), flen[1]/(yglob-1), flen[2]/(zglob-1))
                str = str + "Estimated mem. required for parallel solve = %.3f MB/proc.\n" % nsmem
                ntot = nsmall[0]*nsmall[1]*nsmall[2]

            else:
                str = str + "Fine mesh spacing = %g x %g x %g A\n" % (flen[0]/(n[0]-1), flen[1]/(n[1]-1), flen[2]/(n[2]-1))
                str = str + "Estimated mem. required for sequential solve = %.3f MB\n" % gmem
                ntot = n[0]*n[1]*n[2]

            str = str + "Number of focusing operations = %i\n" % nfocus

            str = str + "\n"
            str = str + "################# ESTIMATED REQUIREMENTS ####################\n"
            str = str + "Memory per processor                   = %.3f MB\n" % (200.0*ntot/1024/1024)
            str = str + "Grid storage requirements (ASCII)      = %.3f MB\n" % (8.0*12*np[0]*np[1]*np[2]*ntot/1024/1024)
            str = str + "Grid storage requirements (XDR)        = %.3f MB\n" % (8.0*ntot*np[0]*np[1]*np[2]/1024/1024)
            str = str + "Time to solve on 667 MHz EV67 Alpha    = %.3f sec\n" % tsolve_alpha
            str = str + "Time to solve on 500 MHz PIII Xeon     = %.3f sec\n" % tsolve_xeon
            str = str + "Time to solve on 400 MHz UltraSparc II = %.3f sec\n" % tsolve_sparc
            str = str + "\n"

        else:
            str = str + "No ATOM entires in file!\n\n"

        return str

def usage():
    psize = Psize()
    usage = "\n"
    usage = usage + "Psize script\n"
    usage = usage + "Usage: psize.py [opts] <filename>\n"
    usage = usage + "Optional Arguments:\n"
    usage = usage + "  --help               : Display this text\n"
    usage = usage + "  --CFAC=<value>       : Factor by which to expand mol dims to\n"
    usage = usage + "                         get coarse grid dims\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("CFAC")
    usage = usage + "  --FADD=<value>       : Amount to add to mol dims to get fine\n"
    usage = usage + "                         grid dims\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("FADD")
    usage = usage + "  --SPACE=<value>      : Desired fine mesh resolution\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("SPACE")
    usage = usage + "  --GMEMFAC=<value>    : Number of bytes per grid point required\n"
    usage = usage + "                         for sequential MG calculation\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("GMEMFAC")
    usage = usage + "  --GMEMCEIL=<value>   : Max MB allowed for sequential MG\n"
    usage = usage + "                         calculation.  Adjust this to force the\n"
    usage = usage + "                         script to perform faster calculations (which\n"
    usage = usage + "                         require more parallelism).\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("GMEMCEIL")
    usage = usage + "  --OFAC=<value>       : Overlap factor between mesh partitions\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("OFAC")
    usage = usage + "  --REDFAC=<value>     : The maximum factor by which a domain\n"
    usage = usage + "                         dimension can be reduced during focusing\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("REDFAC")
    usage = usage + "  --TFAC_ALPHA=<value> : Number of sec/unknown for setup/solve on 667\n"
    usage = usage + "                         MHz EV67 Alpha CPU -- VERY ROUGH ESTIMATE\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("TFAC_ALPHA")
    usage = usage + "  --TFAC_XEON=<value>  : Number of sec/unknown for setup/solve on 500\n"
    usage = usage + "                         MHz PIII Xeon CPU -- VERY ROUGH ESTIMATE\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("TFAC_XEON")
    usage = usage + "  --TFAC_SPARC=<value> : Number of sec/unknown for setup/solve on 400\n"
    usage = usage + "                         MHz UltraSPARC II CPU -- VERY ROUGH ESTIMATE\n"
    usage = usage + "                         [default = %g]\n" % psize.getConstant("TFAC_SPARC")

    
    sys.stderr.write(usage)
    sys.exit(2)

def main():
    import getopt
    filename = ""
    shortOptList = ""
    longOptList = ["help", "CFAC=", "FADD=", "SPACE=", "GMEMFAC=", "GMEMCEIL=", "OFAC=", "REDFAC=", "TFAC_ALPHA=", "TFAC_XEON=", "TFAC_ALPHA="]
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptList, longOptList)
    except getopt.GetoptError, details:
        sys.stderr.write("Option error (%s)!\n" % details)
        usage()
    if len(args) != 1: 
        sys.stderr.write("Invalid argument list!\n")
        usage()
    else:
        filename = args[0]

    psize = Psize()    

    for o, a in opts:
        if o == "--help":
            usage()
        if o == "--CFAC":
            psize.setConstant("CFAC", float(a))
        if o == "--FADD":
            psize.setConstant("FADD", int(a))
        if o == "--SPACE":
            psize.setConstant("SPACE", float(a))
        if o == "--GMEMFAC":
            psize.setConstant("GMEMFAC", int(a))
        if o == "--GMEMCEIL":
            psize.setConstant("GMEMCEIL",  int(a))
        if o == "--OFAC":
            psize.setConstant("OFAC", float(a))
        if o == "--REDFAC":
            psize.setConstant("REDFAC", float(a))
        if o == "--TFAC_ALPHA":
            psize.setConstant("TFAC_ALPHA", float(a))
        if o == "--TFAC_XEON":
            psize.setConstant("TFAC_XEON",  float(a))
        if o == "--TFAC_SPARC":    
            psize.setConstant("TFAC_SPARC",  float(a))

    psize.runPsize(filename)
    sys.stdout.write("Default constants used (./psize.py --help for more information):\n")
    sys.stdout.write("%s\n" % psize.constants)
    sys.stdout.write(psize.printResults())

    
if __name__ == "__main__": main()
