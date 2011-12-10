############################################################################
#
# Authors: Stefano Forli, Ruth Huey
#
# Copyright: M. Sanner TSRI 2010
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/waterMapBuilder.py,v 1.2 2010/09/30 17:13:36 rhuey Exp $
#
# $Id: waterMapBuilder.py,v 1.2 2010/09/30 17:13:36 rhuey Exp $
#

"""
This module implements the WaterMapBuilder class which combines an AutoGrid oxygen map ("receptor.OA.map") 
and an AutoGrid hydrogen affinity map ("receptor.HD.map") to create a water affinity map ("receptor.W.map").

Options and default values include:
default_weight   0.5
mode            "BEST"
entropy          0.0
O_weight         1.0
hd_boost         1.0

"""

class WaterMapBuilder:
    """Class for creating an AutoGrid water map 'stem.W.map' by combining an oxygen map 'stem.OA.map' and a
    hydrogen map 'stem.HD.map'. 
    The available methods for combining the maps based on each pair of values from OA.map and  HD.map include: 
        'best': use the smaller of [OA.map value *self.O_weight, HD.map value *self.hd_boost]
        'avg': if not positive, if either is > use self.ENTROPY =0.0 by default 
               else (OA.map value + HD.map value)/2 ) * self.weight
        'coop': if not positive, if either is > use self.ENTROPY =0.0 by default 
               else ((OA.map value  + HD.map value)) * self.weight
        'boost': if not positive, if either is > use self.ENTROPY =0.0 by default 
                else ((OA.map value *7 + HD.map value *7)/2) * self.weight
        'best_w': if not positive, if either is > use self.ENTROPY =0.0 by default 
                else if OA.map.value < HD.map value :
                    self.O += 1
                    return OA.map.value*self.weight
                else:
                    self.H += 1
                    return HD.map.value*self.weight

    """

    def __init__(self, receptor_file=None, OAmap_file=None, HDmap_file=None, 
                    Wmap_file=None, default_weight=0.5, mix_method="best",
                    ENTROPY=0.0, O_weight=1, hd_boost=1, weight=None, verbose=False):
        #NB: mix_method==mode
        self.receptor_file = receptor_file
        self.OAmap_file = OAmap_file   #firstMap
        self.HDmap_file = HDmap_file   #secondMap
        self.Wmap_file = Wmap_file
        self.default_weight = default_weight
        self.mix_method = mix_method
        self.ENTROPY = ENTROPY
        self.O_weight = O_weight
        self.hd_boost = hd_boost
        self.weight = weight
        if weight is None: 
            weight = default_weight
            self.weight = default_weight
        self.verbose = verbose
        output_stem = None
        if receptor_file is None and (OAmap_file is None or HDmap_file is None):
            raise "receptor filename or both oxygen map file AND hydrogen map file must be specified"
        #???
        self.O = 0
        self.H = 0
        if receptor_file is not None:
            output_stem = os.path.splitext(os.path.basename(receptor_file))[0]
            if self.verbose: print "set output_stem from receptor_file to ", output_stem
            self.OAmap_file = output_stem + ".OA.map"
            self.HDmap_file = output_stem + ".HD.map"
        elif OAmap_file is not None and HDmap_file is not None:
            output_stem = OAmap_file.split('.')[0] 
            if self.verbose: print "set output_stem from OAmap to ", output_stem
        if output_stem is None:
            raise "output_stem is None!"
        self.output_stem = output_stem
        if self.Wmap_file is None:
            self.Wmap_file = "protein.W.map.%s.w%1.2f.O%1.1f.H%1.1f.E%1.1f" % (mix_method, weight, O_weight, hd_boost, ENTROPY )

            #by default self.Wmap_file = output_stem + ".W.map"
            if self.verbose: print "set self.Wmap_file to ", self.Wmap_file
        if self.OAmap_file is None:
            raise "oxygen map file not specified"
        if self.HDmap_file is None:
            raise "hydrogen map file not specified"
        if mix_method not in ['best','avg','coop','boost','best_w']: #mode
            msg = 'Unrecognized mix_method:' + mix_method 
            raise  msg
        if mix_method == 'best':
            self.method = self.best
        elif mix_method == 'avg':
            self.method = self.avg
        elif mix_method == 'coop':
            self.method = self.coop
        elif mix_method == 'boost':
            self.default_weight /= 0.5
            self.method = self.boost
        elif mix_method == 'best_w':
            self.method = self.best_w


    def best(self, first, second, positive=False ): #, hd_boost = hd_boost):
        #global O, H   #@@@@@@@@@@@@@
        """
        - Return the best value between two map values
        - If any positive values are found, 0 is returned by default
        - by default it should be first=OA, second=HD
        """
        #print "in best" 
        first *= self.O_weight
        second *= self.hd_boost
        #print "first, second:", first, second

        if not positive and (first > 0 or second > 0):
                # TODO add a smoothing function for positive values?
                return self.ENTROPY
        else:
            if first < second :
                self.O += 1
                return first * self.weight
            else:
                self.H += 1
                #print "First: %2.6f\tSecond: %2.6f" % (first, second)
                return second * self.weight


    def best_w(self, first, second, positive=False):
        #global O, H   #@@@@@@@@@@@@@
        # TODO test different parameters to check if there's
        # a real advantage
        """
        - Return the best value between two map values
        - If any positive values are found, 0 is returned by default
        - WEIGHTED: the OA map is scaled by weight to account for the higher
                    intensity of OA maps versus HD
        """
        
        first *= self.O_weight
        second *= self.hd_boost

        if not positive and (first > 0 or second > 0):
                return 0 + self.ENTROPY
        else:
            if first < second :
                self.O += 1
                return first*self.weight
            else:
                self.H += 1
                return second*self.weight



    def avg(self, first, second, positive = False):
        """
        - Return the average value between two map values
        - If any positive values are found, 0 is returned by default
        """
        if not positive and (first > 0 or second > 0):
                return self.ENTROPY
        else:
            value = ( (first + second)/2 ) * self.weight
            return value


    def coop(self, first, second, positive = False):
        """
        - Return the sum of value between two map values
        - If any positive values are found, 0 is returned by default
        """
        if not positive and (first > 0 or second > 0):
            return self.ENTROPY
        else:
            value = ((first + second)) * self.weight
            return value


    def boost(self, first, second, positive = False):
        if not positive and (first > 0 or second > 0):
            return self.ENTROPY
        else:
            value =  ((first*7+second*7)/2) * self.weight
            return value


    def build(self):
        if self.verbose:
            print "in build:",
            print " self.OAmap_file=", self.OAmap_file
            print " self.HDmap_file=", self.HDmap_file
            print " self.Wmap_file=", self.Wmap_file
        try:
            #MAP1 = open(firstMap, 'r')
            oa_ptr = open(self.OAmap_file)
            oa_lines = oa_ptr.readlines()
            oa_ptr.close()
            if self.verbose: print "read %d lines from %s"%(len(oa_lines),self.OAmap_file)
        except:
            print "ERROR: unable to open map file: %s" % self.OAmap_file
            exit(1)

        try:
            #MAP2 = open(secondMap, 'r' )
            hd_ptr = open(self.HDmap_file)
            hd_lines = hd_ptr.readlines()
            hd_ptr.close()
            if self.verbose: print "read %d lines from %s"%(len(hd_lines),self.HDmap_file)
        except:
            #print "ERROR: the map %s can't be open" % self.HDmap_file
            print "ERROR: unable to open map file: %s" % self.HDmap_file
            exit(1)

        #sanity checks:
        assert len(oa_lines)==len(hd_lines)

        try:
            if self.verbose: print "opening ", self.Wmap_file
            optr = open(self.Wmap_file, 'w' )
            if self.verbose: print "opened water mapfile: ", self.Wmap_file
            #DIFFERENCE = open(self.Wmap_file, 'w' )
        except:
            print "ERROR: impossible to open the output map", self.Wmap_file
            exit(1)

        #get the header lines from the oxygen map
        #HEADER = lines[:6]
        for ll in oa_lines[:6]:
            optr.write(ll)
            #if self.verbose: print "wrote ", ll
        # process the rest of the lines
        mode = self.mix_method
        if self.verbose: print "using mode =", mode, 'to create W map_file=', self.Wmap_file
        for i in range(6, len(hd_lines)):
            oxy_val=oa_lines[i].strip()
            oxy_val = float(oxy_val)
            hyd_val=hd_lines[i].strip()
            hyd_val = float(hyd_val)
            optr.write("%1.4f\n"%apply(self.method, (oxy_val, hyd_val)))
        optr.close()     
        if self.verbose: print "closed ", self.Wmap_file

#        # output the results
#        for index in range(len(HEADER)):
#            DIFFERENCE.write(str(HEADER[index]))
#####
#        for index in range(len(total)):
#            try:
#                value = "%1.4f" % (total[index])
#            except:
#                print "problem"
#            DIFFERENCE.write(value+'\n')
#        DIFFERENCE.close
        if mode in ['best', 'best_w']:
            O = float(self.O)
            H = float(self.H)
            tot = H+O
            if tot == 0: raise 'O+H == zero !'
            o_pc = (O/tot)*100
            h_pc = (H/tot)*100

            print "\n     output results   "
            print "  --------------------"
            print "  OA points : %3.2f%s" % ( o_pc , "%")
            print "  HD points : %3.2f%s\n" % ( h_pc, "%")

