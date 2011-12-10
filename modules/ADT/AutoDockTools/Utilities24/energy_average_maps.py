#!/usr/bin/env python
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/energy_average_maps.py,v 1.4 2009/03/23 16:40:44 rhuey Exp $
import glob, math

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: energy_average_maps.py "
        print 
        print "    Invoke this script in directory containing all maps..."
        print "    Optional parameters:"
        print "        [-s]    map stem ('weighted')"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 's:vh')
    except getopt.GetoptError, msg:
        print 'energy_average_maps.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    # initialize optional parameters
    #-s: map stem
    map_stem = 'weighted'
    verbose = False
    #'s:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-s', '--s'):
            map_stem = a
            if verbose: print 'set map_stem to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    rt = 0.592
    stems = {}
    atomtypes = {}

    all_maps = glob.glob("*.map")
    for m in all_maps:
        ll = m.split('.')
        stems[ll[0]] = 1
        atomtypes[ll[1]] = 1

    #remove 'e' and 'd' which are computed differently
    del(atomtypes['e'])
    del(atomtypes['d'])

    all_stems = stems.keys()
    all_stems.sort()
    if verbose: print 'all_stems=', all_stems

    all_types = atomtypes.keys()
    all_types.sort()
    if verbose: print 'all_types=', all_types

    all_lines = {}
    for t in all_types:      #loop over atom types
        for s in all_stems:  #loop over receptors
            all_lines[s] = {}
            fn = "%s.%s.map" %(s,t)
            if verbose: print "opening ", fn
            fptr = open(fn)
            all = fptr.readlines()
            fptr.close()
            #skip the first 6 lines of header and convert readlines output(strings) to floats
            all_lines[s][t] = map(float, all[6:])
        num_pts = len(all_lines[s][t]) 
        all_wt_values = []
        for i in range(num_pts):
            vals = []
            wvals = []
            wtot = 0
            eboltz = 0
            for s in all_stems:
                val = all_lines[s][t][i]
                vals.append(val)
                wval = math.e**(-1.*val/rt)
                wvals.append(wval)
                wtot+= wval #add this value
            for j in range(len(all_stems)):
                if wtot!=0:    #if it is, eboltz stays 0
                    eboltz += vals[j] * wvals[j]/wtot
            all_wt_values.append(eboltz)
        #default map_stem is 'weighted'
        fn = map_stem + '.'+ t +'.map'
        fptr = open(fn, 'w')
        for l in all[:6]:
            fptr.write(l)
        for v in all_wt_values:
            l = '%.3f\n'%v
            fptr.write(l)
        fptr.close()
        if verbose: print "wrote ", fn
    #------------------------------------------------------------------
    #   d,e maps
    #------------------------------------------------------------------
    for t in ['e','d']:
        for s in all_stems:  #loop over receptors
            all_lines[s] = {}
            fn = "%s.%s.map" %(s,t)
            if verbose: print "opening ", fn
            fptr = open(fn)
            all = fptr.readlines()
            fptr.close()
            #skip the first 6 lines of header and convert to float
            all_lines[s][t] = map(float, all[6:])
        #perhaps check this is the same as the first part
        #assert len(all_lines[s][t])==num_pts 
        num_pts = len(all_lines[s][t]) 
        all_wt_values = []
        ct = len(all_stems)
        for i in range(num_pts):
            #get list of all the values for this point
            pt_vals = []
            for s in all_stems:
                val = all_lines[s][t][i]
                pt_vals.append(val)
            #get list of absolute value of all the values for this point
            abs_values = map(abs, pt_vals)
            #get mininum of list of absolute values
            mval = min(abs_values)
            #get index of the mininum of list of absolute values
            ind = abs_values.index(mval)
            #get corresponding signed value
            all_wt_values.append(pt_vals[ind])
            #all_wt_values.append(eboltz)
        fn = map_stem + '.'+t+'.map'
        fptr = open(fn, 'w')
        for l in all[:6]:
            fptr.write(l)
        for v in all_wt_values:
            l = '%.3f\n'%v
            fptr.write(l)
        fptr.close()
        if verbose: print "wrote ", fn

    #------------------------------------------------------------------
    #   .fld and .xyz files: use last stem 's'
    #------------------------------------------------------------------
    #xyz file
    xyz_fn = s + '.maps.xyz'
    fptr = open(xyz_fn)
    lines = fptr.readlines()
    fptr.close()
    if verbose: print "read ", xyz_fn
    new_fn = map_stem + '.maps.xyz'
    optr = open(new_fn, 'w')
    for l in lines: 
        optr.write(l)
    optr.close()    
    if verbose: print "wrote ", new_fn
    #fld file
    fld_fn = s + '.maps.fld'
    fptr = open(fld_fn)
    lines = fptr.readlines()
    fptr.close()
    if verbose: print "read ", fld_fn
    new_fn = map_stem + '.maps.fld'
    optr = open(new_fn, 'w')
    for l in lines: 
        optr.write(l.replace(s, map_stem))
    optr.close()    
    if verbose: print "wrote ", new_fn
