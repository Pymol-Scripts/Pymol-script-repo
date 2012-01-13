#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/trajParser.py,v 1.2 2006/10/17 18:22:13 annao Exp $
#
# $Id: trajParser.py,v 1.2 2006/10/17 18:22:13 annao Exp $
#

## This module implements classes for parsing GROMACS trajectory files

from struct import calcsize
from xdrlib import Unpacker
DIM = 3
import os

class trrParser:
    """ Parses .trr Gromacs trajectory file """
    
    def __init__(self, file):

        self.nframes = 0
        self.headers=[]
##         header = {
##             "ir_size":None,
##             "e_size": None,
##             "box_size" :None,
##             "vir_size": None,
##             "pres_size": None,
##             "top_size ": None,
##             "sym_size ": None,
##             "x_size ": None,
##             "v_size ": None,
##             "f_size ": None,
##             "natoms ": None,
##             "step ": None,
##             "nre ": None,
##             "version":'GMX_trn_file',
##             "magicnum" :1993,
##             "bDouble ": None,
##             }
        self.file = file
        self.coords = []
        self.velocities = {}
        self.forces = {}


    def nFloatSize(self, h):

        nflsize=0;
        if h["box_size"]:
            nflsize = h["box_size"]/(DIM*DIM);
        elif h["x_size"]:
            nflsize = h["x_size"]/(h["natoms"]*DIM);
        elif h["v_size"]:
            nflsize = h["v_size"]/(h["natoms"]*DIM);
        elif h["f_size"]:
            nflsize = h["f_size"]/(h["natoms"]*DIM);
        else: 
            print "Can not determine precision of trr file"
  
        if (nflsize != calcsize("f")) and (nflsize != calcsize("d")):
            print "Float size %d. Maybe different CPU?"%nflsize
      
        return nflsize

    def read(self):
        fext = os.path.splitext(self.file)[-1]
        assert fext == ".trr"
        fp = open(self.file, "rb")
        self.data = data = fp.read()
        self.coords = []
        self.v = {}
        self.f = {}
        self.up = Unpacker(data)
        curpos = self.up.get_position()
        datasize = len(data)
        nframe = 0
        #each frame begins with a header
        while curpos < datasize:
            #print "current position:", curpos
            h = self.readHeader(nframe)
            self.headers.append(h)
            self.readData(nframe)
            nframe = nframe + 1
            curpos = self.up.get_position()
        #print "end of readTraj, cur position : %d, datazize: %d" %(self.up.get_position(), datasize)
        self.nframes = nframe
        if self.nframes:
            return 1
        else:
            return 0

    def readHeader(self, nframe):
        #print "reading header, frame %d" %nframe
        up = self.up
        header = {}
        header["magicnum"] = up.unpack_int()
        #print "magicnum: ", header["magicnum"]
        #1993
        i1 = up.unpack_int()
        #13
        version = up.unpack_string()
        #'GMX_trn_file'
        header["version"] = version
        #print "version:", version
        header["ir_size"] = up.unpack_int()
        #print "ir_size=", header["ir_size"]
        header["e_size"]=up.unpack_int()
        #print "e_size=", header["e_size"]
        header["box_size"] =up.unpack_int()
        #print "box_size=", header["box_size"]
        header["vir_size"] = up.unpack_int()
        #print "vir_size=", header["vir_size"]
        header["pres_size"] = up.unpack_int()
        #print "pres_size=", header["pres_size"]
        header["top_size"]=up.unpack_int()
        #print "top_size=", header["top_size"]

        header["sym_size"]=up.unpack_int()
        #print "sym_size=", header["sym_size"]

        header["x_size"]=up.unpack_int()
        #print "x_size=", header["x_size"]

        header["v_size"]=up.unpack_int()
        #print "v_size=", header["v_size"]

        header["f_size"]=up.unpack_int()
        #print "f_size=", header["f_size"]

        header["natoms"]=up.unpack_int()
        #print "natoms=", header["natoms"]
        header["step"]=up.unpack_int()
        #print "step=", header["step"]
        header["nre"]=up.unpack_int()
        #print "nre=", header["nre"]

        
        if self.nFloatSize(header) == calcsize("d"):
            header["bDouble"] = True
        else:
            header["bDouble"] = False

        if header["bDouble"]:
            header["time"] = up.unpack_double()
            header["lam"] = up.umpack_double()
        else:
            header["time"] = up.unpack_float()
            header["lam"] = up.unpack_float()

        #print "time = ", header["time"]
        #print "lambda = ", header["lam"]
        #print "natoms=%10d  step=%10d  time=%10g  lambda=%10g"% (header["natoms"],header["step"],header["time"],header["lam"])
        #print "current position:",  up.get_position()
        return header


    def readData(self, nframe):
        up = self.up
        h = self.headers[nframe]
        box = []
        if h["box_size"] != 0 :
            for i in range(3):
                box.append(up.unpack_farray(3, up.unpack_float))
            #print " box (3x3):"
            #print box
            self.headers[nframe]["box"] = box
        pv = []
        if h["vir_size"] != 0:
            for i in range(3):
                pv.append(up.unpack_farray(3, up.unpack_float))
            #print "pv:"
            #print pv
            self.headers[nframe]["pv"] = pv
        if h["pres_size"]!= 0:
            pv.append(up.unpack_farray(3, up.unpack_float))
            #print "pv:"
            #print pv
            self.headers[nframe]["pv"] = pv
        natoms = h["natoms"]
        if h["x_size"] != 0:
            x= []
            for i in range (natoms): 
                x.append(up.unpack_farray(3, up.unpack_float))
            #self.coords["frame%d"%nframe] = x
            self.coords.append(x)

        if h["v_size"] != 0:
            v = []
            for i in range (natoms):
                v.append(up.unpack_farray(3, up.unpack_float))
            self.velocities["frame%d"%nframe] = v
        if h["f_size"] != 0:
            f = []
            for i in range (natoms):
                f.append(up.unpack_farray(3, up.unpack_float))
            self.forces["frame%d"%nframe] = f
        #print up.get_position()





class xtcParser:

    """ Parses .xtc Gromacs trajectory file """

    def __init__(self, file):
        
        self.nframes = 0
        self.headers=[]
##         header = {
##             'step ': None,
##             'frame':None,
##             'time':None,
##             'prec':None
##             }
        self.file = file
        self.coords = []
        self.headers = []
        self.nframes = 0
        self.velocities = None #not available inthis file format
        self.forces = None
        # try to import a platform-dependent module xtcparcer
        # containing C functions to read this file format
        try:
            from cMolKit import xtcparser
        except:
            print "WARNING: could not import cMolKit.xtcparser - No parser is available for xtc files."
            self.file=None


    def read(self):
        if not self.file:
            return 0
        assert os.path.exists(self.file)
        fext = os.path.splitext(self.file)[-1]
        assert fext == ".xtc"
        
        from cMolKit import xtcparser
        self.coords, self.headers = xtcparser.read_xtc(self.file)
        if self.coords:
            self.nframes = len(self.coords)
            return 1
        else:
            return 0
