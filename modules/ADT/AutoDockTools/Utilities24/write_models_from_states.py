#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_models_from_states.py,v 1.4 2010/04/05 21:07:42 rhuey Exp $
#
# $Id: write_models_from_states.py,v 1.4 2010/04/05 21:07:42 rhuey Exp $
#
import os, sys, sys
from string import join
from MolKit import Read
from mglutil.math.statetocoords import StateToCoords
from AutoDockTools.Conformation import Conformation

from AutoDockTools.Docking import Docking


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: write_modes_from_states.py -l ligandfile -s statefile -o multimodelfile"
        print
        print "    Description of command..."
        print "         -l     ligandfile name"
        print "         -s     statefile name"
        print "    Optional parameters:"
        print "        [-o]    multimodel output filename "
        print "        [-S]    single string replacing statefile eg:"
        print "    'State: 29.303 14.415 23.603 0.5609 0.4518 0.2662 -0.6406 -20.89 -0.65 81.86 -17.36 28.83 -10.80 -23.98 114.21'"
        print "        [-e]    statefile includes energy"
        print "        [-z]    use zero origin"
        print "        [-i]    interim state->apply quaternion before 'about' translation"
        print "        [-n]    index of energy on energy line: default is 8"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:S:s:o:n:ezivh')
    except getopt.GetoptError, msg:
        print 'write_modes_from_states.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: ligandfile name
    ligandfile =  None
    #-S: state
    SINGLESTATE =  None #0? or False?
    #-s: statefile name
    statefile =  None
    #-o multimodel outputfilename
    outputfile = None
    #-e states_have_energy
    states_have_energy = 0
    #-n index of energy
    index_of_energy = 8
    #-z use_zero_origin 
    use_zero_origin = False
    #-i interim_state 
    interim_state = False


    # initialize optional parameter
    #-v verbose best only
    verbose = False

    #'s:o:evh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-l', '--l'):
            ligandfile = a
            if verbose: print 'set ligandfile to ', a
        if o in ('-S', '--S'):
            SINGLESTATE = 1
            if verbose: print 'set SINGLESTATE to ', a
            #['30.691_15.206_23.914_-0.3807_0.0201_0.9215_0.0747_-7.06_93.27_127.02_-130.67_7.04_-95.44_-4.91_-126.85']"
            if args[0].find('_')>-1:
                args[0] = args[0].replace("_"," ")
            states = ["State: " + join(args)]
        if o in ('-s', '--s'):
            statefile = a
            if verbose: print 'set statefile to ', a
        if o in ('-o', '--o'):
            outputfile = a
            if verbose: print 'set outputfile to ', a
        if o in ('-e', '--e'):
            states_have_energy = True
            if verbose: print 'set states_have_energy to ', states_have_energy
        if o in ('-n', '--n'):
            index_of_energy = int(a)
            if verbose: print 'set index_of_energy to ', index_of_energy
        if o in ('-z', '--z'):
            use_zero_origin = True
            if verbose: print 'set use_zero_origin to ', use_zero_origin
        if o in ('-i', '--i'):
            interim_state = True
            if verbose: print 'set interim_state to ', interim_state
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not ligandfile:
        print 'write_modes_from_states.py: ligandfile must be specified.'
        usage()
        sys.exit()

    if not statefile and not SINGLESTATE:
        print 'write_modes_from_states.py: SINGLESTATE or statefile must be specified.'
        usage()
        sys.exit()

    if not outputfile:
        if verbose: print 'write_modes_from_states.py: outputfile not specified. Using stdout'
        #usage()
        #sys.exit()


    lig = Read(ligandfile)
    if not len(lig):
        print "no ligand found in ", ligandfile
        sys.exit()
    lig = lig[0]
    if not hasattr(lig, 'ndihe'):
        print ligandfile + "molecule has no torsion tree"
        sys.exit()
    lig.buildBondsByDistance()
    # add extra slot to ._coords for changing coordinates
    lig.allAtoms.addConformation(lig.allAtoms.coords)
    #?is this necessary
    lig.allAtoms.setConformation(1)
    ntors = lig.ndihe
    length_of_state = 7+lig.ndihe
    # @@ handle to the input ligLines
    ligLines = lig.parser.allLines

    #setup StateToCoords object
    origin = lig.getCenter()
    if use_zero_origin or interim_state:
        origin = [0.,0.,0.]    
    #note: index of _coords to use is always 1
    lig.stoc = StateToCoords(lig, origin, 1)
    outptr = sys.stdout
    if outputfile:
        outptr = open(outputfile, 'w')


    #if SINGLESTATE:
        # eg:
        #"State: 29.303 14.415 23.603 0.5609 0.4518 0.2662 -0.6406 -20.89 -0.65 81.86 -17.36 28.83 -10.80 -23.98 114.21"
        #states = state.split()
        #['State:', '29.303', '14.415', '23.603', '0.5609', '0.4518', '0.2662', '-0.6406', '-20.89', '-0.65', '81.86', '-17.36', '28.83', '-10.80', '-23.98', '114.21']

    if statefile:
        sptr = open(statefile)
        states = sptr.readlines()
        sptr.close()
        if not len(states):
            print "no states found in ", statefile
            sys.exit()

    state_list = []
    ctr = 1
    count = len(states)
    if not SINGLESTATE and states_have_energy:
        for i in range(0,len(states),2):
            sline = states[i]
            eline =states[i+1]
            #build a states from each of the lines in statefile
            #State:\t  4.847  -2.386  14.760  -0.413  0.552 -0.724  4.257     58.27  -33.47  -87.92  134.64  -36.46  114.79  -44.86  -74.96 -118.53   77.29  139.08   78.23  -52.09  -12.69   35.08 -118.21 -175.94\n'
            fl = map(float, sline.split()[1:])
            assert len(fl)==length_of_state
            # 0  1  2  3  4  5  6 [7....
            #[t1,t2,t3,q1,q2,q3,q4,tors1,tors2, tors3....
            #
            translation = fl[:3]
            quaternion = [fl[6],fl[3:6]]
            #energy = eline.split()[8]
            energy = eline.split()[index_of_energy]
            torsion_angles = []
            if ntors>0:
                torsion_angles = fl[7:]
                newConf = Conformation(lig,origin,translation, quaternion, torsion_angles)
                newCrds = newConf.getCoords()
                if interim_state:
                    #here's where to add back the origin or ?
                    newCrds -= origin
                #write some MODEL stuff then newCrds
                ostr = "MODEL %d\n"%ctr
                outptr.write(ostr)
                ctr += 1
                ostr = "REMARK AD4 RESULT: %s\n" %energy #put energy here...
                outptr.write(ostr)
                #lig.parser.write_with_new_coords(newCrds,outptr)
                ct = 0
                for l in ligLines:
                    if l.find("ATOM")==0 or l.find("HETATM")==0:
                        cc = newCrds[ct]
                        ct = ct + 1
                        new_l = l[:30]+"%8.3f%8.3f%8.3f" %(cc[0],cc[1],cc[2]) + l[54:]
                    else:
                        new_l = l
                    outptr.write(new_l) 
                if verbose: print "wrote ", outputfile
                ostr = "ENDMDL\n"
                outptr.write(ostr)
                i+=1
    else:
        for sline in states:
            #build a state from each of the lines in states
            #State:\t  4.847  -2.386  14.760  -0.413  0.552 -0.724  4.257     58.27  -33.47  -87.92  134.64  -36.46  114.79  -44.86  -74.96 -118.53   77.29  139.08   78.23  -52.09  -12.69   35.08 -118.21 -175.94\n'
            fl = map(float, sline.split()[1:])
            assert len(fl)==length_of_state
            # 0  1  2  3  4  5  6 [7....
            #[t1,t2,t3,q1,q2,q3,q4,tors1,tors2, tors3....
            #
            translation = fl[:3]
            # interpret axis_angle as quaternion x y z w
            # use quaternion as w, (x,y,z) for mglutil/math/transformation.py class
            quaternion = [fl[6],fl[3:6]]
            torsion_angles = []
            if ntors>0:
                torsion_angles = fl[7:]
                newConf = Conformation(lig,origin,translation, quaternion, torsion_angles)
                newCrds = newConf.getCoords()
                #write some MODEL stuff then newCrds
                ostr = "MODEL %d\n"%ctr
                outptr.write(ostr)
                ctr += 1
                ostr = "REMARK AD4 RESULT: n/a\n" #put energy here...
                outptr.write(ostr)
                #lig.parser.write_with_new_coords(newCrds,outptr)
                ct = 0
                for l in ligLines:
                    if l.find("ATOM")==0 or l.find("HETATM")==0:
                        cc = newCrds[ct]
                        ct = ct + 1
                        new_l = l[:30]+"%8.3f%8.3f%8.3f" %(cc[0],cc[1],cc[2]) + l[54:]
                    else:
                        new_l = l
                    outptr.write(new_l) 
                if verbose: print "wrote ", outputfile
                ostr = "ENDMDL\n"
                outptr.write(ostr)
    if verbose: print "Done!"
    outptr.close()            

# To execute this command type:
# write_modes_from_states.py -d docking_filename 
# optional arguments
# -o outputfile_stem (default is ligandname)
