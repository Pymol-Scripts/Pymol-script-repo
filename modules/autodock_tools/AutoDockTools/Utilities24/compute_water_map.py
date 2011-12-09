#!/usr/bin/env python
#
# Combine OA and HD maps for providing the W map
# for water molecules affinity
# The first map should be OA, the second HD
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_water_map.py,v 1.1 2010/09/30 17:21:36 rhuey Exp $
#
# $Id: compute_water_map.py,v 1.1 2010/09/30 17:21:36 rhuey Exp $
#
# Authors: Stefano Forli, Ruth Huey
#
#

from AutoDockTools.waterMapBuilder import WaterMapBuilder

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print """\n
                                        /\/\   __ _ _ __  
                                       /    \ / _` | '_ \ 
                                      / /\/\ \ (_| | |_) |           .-.
                                      \/    \/\__,_| .__/           (   )_
                                                   |_|              _/-'(_)
                                   __    __      _                 (_)
                                  / / /\ \ \__ _| |_ ___ _ __ 
                                  \ \/  \/ / _` | __/ _ \ '__|
                                   \  /\  / (_| | ||  __/ |   
                                    \/  \/ \__,_|\__\___|_|   
                                  
        """
        print "\tINPUT\n\t\tEither the receptor filename or the map files.\n"
        print "\t\t\t-r protein.pdbqt    :   the receptor filename is used to guess the OA and HD filenames"
        print "\n\t\t\t\t *OR*\n"
        print "\t\t\t-o protein.OA.map   :   OA map filename (omit if -r is used)"
        print "\t\t\t-h protein.HD.map   :   HD map filename (omit if -r is used)"
        print "\tOUTPUT\n\t\tW map.\n"
        print "\tOPTIONS\n"
        print "\t\t-m (best,avg,coop..):   mix method (default 'best')"
        print "\t\t                        (Other methods are 'avg','coop', 'boost','best_w')"
        print "\t\t-e float            :   entropy penalty value (default: 0.0)"
        print "\t\t-O float            :   oxygen map weight (default: 1.0)"
        print "\t\t-H float            :   hydrogen map weight (default: 1.0)"
        print "\t\t-s protein.W.map    :   W map output filename" 
        print "\t\t   (default: protein.W.map.[mix method].w[eight].o[xygen weight].h[ydrogen weight].E[ntropy correction])"
        print "\n\n"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'o:h:m:e:O:H:s:w:r:vU')
        # -o protein.OA.map
        # -h protein.HD.map
        # -m mix_method best,avg,coop...
        # -e entropy
        # -O oxygen_weight
        # -H hydrogen_boost
        # -s protein.W.map
        # -w weight 0.5
        # -r protein.pdbqt
        # -v verbose
        # -U print Usage statement
    except getopt.GetoptError, msg:
        print 'compute_water_map.py: %s' %msg
        usage()
        sys.exit(2)
        #exit(2) #???

    # initialize required parameters
    # -o protein.OA.map
    oa_map = None
    # -h protein.HD.map
    hd_map = None
    # -m mix_method or mode
    mix_method = 'best'  #best,avg, coop, boost
    # -e entropy 
    entropy = 0.0
    # -O oxygen_weight
    O_weight = 1.0
    # -H hd_boost
    hd_boost = 1.0
    # -s protein.W.map
    w_map = None
    # -w weight
    weight = 0.5
    default_weight = weight
    # -r protein.pdbqt
    protein = None

    ENTROPY = 0.0
    verbose = False

    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            protein =  os.path.basename(a)
            if verbose: print "set protein to ",a
            name = os.path.splitext(protein)[0]
            #name = protein.rsplit(".")[0]
            oxygenMap = name+".OA.map"
            hydrogenMap = name+".HD.map"
            if verbose:
                print "      OA map ->", oxygenMap
                print "      HD map ->", hydrogenMap
        if o in ('-o', '--o'):
            oa_map = a
        if o in ('-h', '--h'):
            hd_map = a

        if o in ('-w','--w'):
            try:
                default_weight = float(a)
                weight = default_weight
            except:
                print "exception using default_weight", a
                exit(1)
        if o in ('-m', '--m'):
            mix_method = a
            if verbose: print "mix_method=", mix_method
            if mix_method == "BOOST":
                default_weight /= .05

        if o in ("-e", "--e"):
            ENTROPY = float(a)

        if o in ("-O", "--O"):
            O_weight = float(a)

        if o in ("-H", "--H"):
            hd_boost = float(a)


        if o in ("-s","--s"):
            #outputFile = a
            w_map = a
        
        if o in ("-U","--U"):
            usage()
            sys.exit()


    if protein is None and oa_map is None and hd_map is None:        
        usage()
        sys.exit()

    if w_map is None:        
        w_map = "protein.W.map.%s.w%1.2f.O%1.1f.H%1.1f.E%1.1f" % (mix_method, weight, O_weight, hd_boost, ENTROPY )
        print "using default outputfile name: ", w_map
        # protein.W.map.BEST.w0.5.O0.3.H1.3.E0


    # initialize required parameters
    print "\n  MapWater calculator\n ====================="
    print "  mode      :  ", mix_method
    print "  weight    :  ", default_weight
    print "  OA_weight :  ", O_weight
    print "  hd_boost  :  ", hd_boost
    print "  Entropy   :  ", ENTROPY
    print "  output    :  ", w_map


    # optional parameters
    #    if o in ('-h', '--'):
    #        usage()
    #        sys.exit()

    if weight==default_weight and verbose:
        print " => Water map weight : DEFAULT [ %1.2f ]" % weight
    wmb = WaterMapBuilder(OAmap_file=oa_map, HDmap_file=hd_map, Wmap_file=w_map,default_weight=default_weight, weight=weight, mix_method=mix_method)
    wmb.build()
# To execute this command type:
# compute_water_map.py 
