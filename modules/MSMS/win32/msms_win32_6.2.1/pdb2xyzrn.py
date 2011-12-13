def pdb2xyzrn(self,pdb_fn,xyzr_fn,xyzrn_fn,use_explicit_h_r=False):
    """ convert pdb file to .xyzr and .xyzrn files. 
        The script pdb_to_xyzr and pdb_to_xyzrn work well in Linux, with some bugs eg for MN and HG.
        This function is designed base on them. 
        The comments are taken from the script to make them as consistent as possible.
        @param use_explicit_h_r: this is the -h parameter 
        @param type: boolean 
    """
    if use_explicit_h_r:    h_select = 4
    else:                   h_select = 5 
    fh = open(pdb_fn)
    fd = fh.readlines()
    fh.close()
    # read the radius table
    # assume the atmtypenumbers file is in the work directory of msms
    atnfile = os.path.normpath('%s/atmtypenumbers' % (self.msms_wd,))
    #print 'atmtypenumbers', atnfile
    fh1 = open(atnfile)
    fd1 = fh1.readlines()
    fh1.close()
    explicit_rad = {}
    united_rad = {}
    npats = 0    # number of patterns 
    respat = {}  # residue patterns (for regular expression)
    atmpat = {}  # atom patterns
    atmnum = {}  # atom number

    for line in fd1:

        if line.startswith('#') or len(line.strip()) == 0: 
            continue
        elif line.startswith('radius'):
            buf = line.strip().split()
            n, explicit_rad[n] = int(buf[1]),float(buf[3])
            if len(buf) <= 4 or buf[4]=='#':  united_rad[n] = explicit_rad[n]
            else:                             united_rad[n] = float(buf[4])
        else:
            buf = line.strip().split()
            respat[npats] = buf[0]
            if respat[npats] == '*': respat[npats] = '.*'
            respat[npats] = '^' + respat[npats] + '$' 
            atmpat[npats] = '^' + buf[1] + '$'
            atmpat[npats].replace('_', ' ')
            atmnum[npats] = int(buf[2])
            if atmnum[npats] not in explicit_rad:
                # the key has no radius --- complain and fake one
                print 'pdb_to_xyzr: error in library file', atnfile, \
                      'entry', buf[0], buf[1], buf[2], 'has no corresponding radius value'
                explicit_rad[atmnum[npats]] = 0.01
                united_rad[atmnum[npats]] = 0.01

            npats += 1

    H_pat  = re.compile('[ \d][HD]')
    output_buf_xyzr  = []
    output_buf_xyzrn = []
    line_num = 0
    for line in fd:

        line_num += 1

        if (line.startswith('ATOM  ') or line.startswith('HETATM')) and len(line) > 54:

            aname   = line[12:16]
            # normally, H atoms should be in the format as explained:
            # http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html
            # but sometimes false format is used, like 1HE2
            if (re.search(H_pat,aname[:2]) is not None) or (aname[0] in {'H':1, 'D':2} and len(aname.strip())==4):
                aname='H'
            # an emprical rule for judging HG or HE atoms
            # for hydrogen atoms, the letter H is put at column 13 only if
            # the length of the atom name is 4
            # http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html
            if aname[:2]=='HG' and len(aname.strip()) < 4: aname = 'HG'
            if aname[:2]=='HE' and len(aname.strip()) < 4: aname = 'HE'
            if aname[:2]=='HF' and len(aname.strip()) < 4: aname = 'HF'

            aname = aname.strip().replace(' ','') # there are spaces in some atoms names: 1sct HEM 
            resname = line[17:20].strip()
            resnum  = line[22:26].strip()  # icode is not consider here!
            x,y,z   = line[30:38], line[38:46], line[46:54]

            for pat in xrange(npats):
                if re.search(atmpat[pat],aname) and re.search(respat[pat], resname):
                    if atmnum[pat] != 15: break # break only if it does not match a H atom pat
                                        # otherwise continue the search, 
                                        #as aname might match a heavy metal pattern.

            if pat == npats:   # Not found
                print "pdb_to_xyzr: error, file",pdb_fn,"line",line_num,"residue", resnum,\
                      "atom pattern",resname, aname,"was not found in ", atnfile 
                output_buf_xyzr.append('%s %s %s 0.01\n' % (x,y,z))
                output_buf_xyzrn.append('%s %s %s 0.01\n' % (x,y,z))
            else:
                output_buf_xyzr.append(
                    '%s %s %s %.2f\n' % \
                    (x,y,z,h_select==5 and united_rad[atmnum[pat]] or explicit_rad[atmnum[pat]]))
                output_buf_xyzrn.append(
                    '%.6f %.6f %.6f %.6f %d %s_%s_%d\n' % \
                    (float(x),float(y),float(z),
                     h_select==5 and united_rad[atmnum[pat]] or explicit_rad[atmnum[pat]],
                     1, aname,resname,int(resnum.strip())))    
        
    print 'output_fn_xyzr', xyzr_fn
    fh = open(xyzr_fn,'w')
    fh.writelines(output_buf_xyzr)
    fh.close()
    print 'output_fn_xyzrn',xyzrn_fn
    fh = open(xyzrn_fn,'w')
    fh.writelines(output_buf_xyzrn)
    fh.close()
    return
