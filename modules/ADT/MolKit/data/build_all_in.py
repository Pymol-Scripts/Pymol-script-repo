############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2001
#
#############################################################################

# $Header: /opt/cvs/python/packages/share1.5/MolKit/data/build_all_in.py,v 1.4 2003/08/29 17:50:16 sophiec Exp $
#
# $Id: build_all_in.py,v 1.4 2003/08/29 17:50:16 sophiec Exp $
#


from string import split, strip


def getLoopInfo(i, allLines):
    #returns a list of lists of pairs of atom names which define loop bonds
    done = 0
    loopList = []
    #to start: i is first line with loop info
    while not done:
        ll = split(allLines[i])
        #print 'gli', ll
        if len(strip(allLines[i]))==0:
            #print 'returning ', allLines[i+1], ' from loopInfo'
            done = 1
        else:
            if ll[0]!='LOOP':
                loopList.append(ll)
            i = i + 1
    #NB: return i for non-blank line after blank line
    return i+1, loopList

        
def getChargeInfo(i, allLines):
    done = 0
    newcharges = []
    ll = split(allLines[i])
    #to start: i is first line with charges
    while not done:
        for ch in ll:
            newcharges.append(float(ch))
        i = i + 1
        ll = split(allLines[i])
        if len(strip(allLines[i]))==0:
            done = 1
    #NB: return i for non-blank line after blank line
    return i+1, newcharges


def getImpropInfo(i, allLines):
    #at this point should be in IMPROPER section
    impropTors = []
    ll = split(allLines[i])
    if ll[0]=='IMPROPER':
        i = i+1
        ll = split(allLines[i])

    while len(ll)==4:
        impropTors.append(ll)
        i = i + 1
        ll = split(allLines[i])
    #NB: return i for non-blank line after blank line
    return i+1, impropTors


def getDUMMList(i, allLines):
    DUMMLIST = []
    #newD['DUMM'] = []
    ll = split(allLines[i])
    if ll[1]=='DUMM':
        #if find DUMM assume there will be three of them
        DUMMLIST.append(ll)
        DUMMLIST.append(split(allLines[i+1]))
        DUMMLIST.append(split(allLines[i+2]))
        i = i + 3
    return DUMMLIST, i


def getAtom(ll):
    # called for each of the atoms:
    #I, IGRAPH,ISYMBL, ITREE, NA, NB, NC,R,THETA, PHI, CHG
    atname = ll[1]
    d = {}
    #keep all amber info
    d['I'] = int(ll[0])
    d['type'] = ll[2]
    d['tree'] = ll[3]
    #d['angleInfo'] = ll[4:7]
    #NA(I)      The atom number to which atom I is connected.
    #Read but ignored for internal coordinates; If cartesian
    #coordinates are used, this must be omitted.
    d['NA'] = int(ll[4])
    #NB(I)      The atom number to which atom I makes an angle along
    #with NA(I).
    #Read but ignored for internal coordinates; If cartesian
    #coordinates are used, this must be omitted.
    d['NB'] = int(ll[5])
    #NC(I)      The atom number to which atom I makes a dihedral along
    #with NA(I) and NB(I).
    #Read but ignored for internal coordinates; If cartesian
    #coordinates are used, this must be omitted.
    d['NC'] = int(ll[6])
    #R(I)       If IFIXC .eq. 'CORRECT' then this is the bond length
    #between atoms I and NA(I)
    #If IFIXC .eq. 'CHANGE' then this is the X coordinate
    #of atom I
    #this is the same as 'blen'
    d['blen'] = float(ll[7])
    #THETA(I)   If IFIXC .eq. 'CORRECT' then it is the bond angle
    #between atom NB(I), NA(I) and I
    #If IFIXC .eq. 'CHANGE' then it is the Y coordinate of
    #atom I
    #this is the same as 'angle'

    #PHI(I)     If IFIXC .eq. 'CORRECT' then it is the dihedral angle
    #between NC(I), NB(I), NA(I) and I
    #If IFIXC .eq. 'CHANGE' then it is the Z coordinate of
    #atom I
    #this is the same as 'torsion'

    #CHRG(I)    The partial atomic charge on atom I
    #this is the same as 'charge'
    if len(ll)==11:
        d['angle'] = float(ll[8])
        d['torsion'] = float(ll[9])
        d['charge'] = float(ll[10])
    elif len(ll)==10:
        d['angle'] = float(ll[8])
        d['torsion'] = float(ll[9])
        d['charge'] = None
    return atname, d


def buildRes(i, allLines):
    #allLines[i] at start is residue NAME
    done = 0
    newD = {}
    newD['NAMRES'] = strip(allLines[i])
    #print 'name=', newD['NAMRES']
    
    #skip blank line and get 'NAMRES, INTX, KFORM'
    i = i + 2
    ll = split(allLines[i])
    #eg ['ALA', 'INT', '1']
    name = ll[0]
    #INTX flg for type of coords to be saved 'INT'internal or 'XYZ'cartesian
    #KFORM: format of output of individual residue files:0 
    #   0 for formatted; 1 for binary 
    newD['INTX,KFORM'] = ll[1:]
    # IFIXC, IOMIT, ISYMDU, IPOS
    i = i + 1
    ll = split(allLines[i])
    #IFIXC input geom type is CORRECT/CHANGE
    #IOMIT OMIT/NOMIT dummy atoms
    #ISYMDU symbol for dummy atoms, 'DU' is preferred
    #IPOS ALL/BEG flag for position of dummay atoms to be deleted
    newD['IFIXC,IOMIT,ISYMDU,IPOS']=ll
    ifixc = ll[0]
    #CUT
    i = i + 1
    ll = split(allLines[i])
    #CUT is cutoff distance for loop closing bonds which cannot be defined by
    #tree structure. Any pair of atoms within this distance is assumed to be
    #bonded.  We recomment that CUT be set to 0.0 and explicit loop closing
    #bonds be defined below.
    newD['CUT'] = ll

    #build residue DUMMLIST, 
    i = i + 1
    newD['DUMM'], i = getDUMMList(i, allLines)

    #build dictionary for each atom 
    endAtoms = 0
    atNameList = newD['atNameList'] = []
    while not endAtoms:
        if len(split(strip(allLines[i])))<9:
            endAtoms = 1
            if ifixc == 'CHANGE' or len(split(strip(allLines[i+1])))==0:
                i = i+1
        else:
            l = strip(allLines[i])
            if not len(l):
                i = i+1
            atName, newD[atName] = getAtom(split(strip(allLines[i])))
            atNameList.append(atName)
        i = i + 1

    ll = split(strip(allLines[i]))
    i = i+1
    # IOPR  specifies what optional additional info follows:
    done = 0
    if not len(ll):
        i = i+1
        ll = split(strip(allLines[i]))
    iopr = ll[0]
    hasLoop = 0
    hasImprop = 0
    hasCharges = 0
    cutoff = len(allLines)-1
    while not done and i < cutoff:
        #print 'not done using', iopr
        if iopr=='LOOP':
            #print ll
            i, newD['loopList'] = getLoopInfo(i, allLines)
            #iopr = split(strip(allLines[i]))[0]
            l = strip(allLines[i])
            if len(l):
                iopr = split(l)[0]
            else:
                i = i+1
                iopr = split(strip(allLines[i]))[0]
            #print 'new iopr after loop =', iopr
            hasLoop = 1
        elif iopr=='CHARGE':
            i, newcharges = getChargeInfo(i, allLines)
            #for CYM overwrite the charge field for 10 atoms
            #put these 10 new values into atom
            lenNewC = len(newcharges)
            if lenNewC!=0:
                for q in range(lenNewC):
                    #order of atoms in atNameList is the same 
                    #as order of newcharges
                    newD[atNameList[q]]['charge']=newcharges[q]
            #iopr = split(strip(allLines[i]))[0]
            l = strip(allLines[i])
            if len(l):
                iopr = split(l)[0]
            else:
                i = i+1
                iopr = split(strip(allLines[i]))[0]
            hasCharges = 1
        elif iopr=='IMPROPER':
            #print 'calling improper on', i
            i, newD['impropTors'] = getImpropInfo(i, allLines)
            l = strip(allLines[i])
            if len(l):
                iopr = split(l)[0]
            else:
                i = i+1
                iopr = split(strip(allLines[i]))[0]
            hasImprop = 1
        elif iopr=='DONE': 
            done = 1        
            if not (hasLoop or hasImprop or hasCharges):
                #print 'in hack for CIP/CIM'
                i = i-1
        elif iopr=='STOP':
            done = 1
        if i==cutoff: done=1
    return i+1, name, newD 


def buildResDict( inFile = 'all_amino94.in'):
    #all.in is alternative dict
    f = open(inFile)
    allLines = f.readlines()
    f.close()
    outdict = {}
    outdict['IDBGEN,IREST,ITYPF'] = split(allLines[0])
    outdict['NAMDBF'] = strip(allLines[1])
    i = 2
    done = 0
    while not done:
        i, name, outdict[name] = buildRes(i, allLines)
        if i>=len(allLines): done = 1
        else:
            ll = split(allLines[i])
        if ll[0]=='STOP': done=1
    outdict['filename'] = inFile
    return outdict


def write(dict, outfile='test_dat.py',dname='all_amino94_dat'):
    fptr = open(outfile,'w')
    outstr = dname + ' = {\n'
    fptr.write(outstr)
    #NB each v is also a dictionary
    for k, v in dict.items():
        if type(v)==type(dict):
            fptr.write('\"'+k+'\"' + ': {')
            for key, val in v.items():
                outstr = '\t'+ '\"'+ key +'\"' + ':' + repr(val) + ',\n'
                fptr.write(outstr)
            fptr.write('},\n')
        else:
            outstr = '\"'+ k + '\"' + ':' + repr(v) + ',\n'
            fptr.write(outstr)
    fptr.write('}')
    fptr.close()


dict=unint_dat = buildResDict('unint.in')
write(dict, outfile='unint_dat.py',dname='unint_dat')
dict=opls_unint_dat = buildResDict('opls_unint.in')
write(dict, outfile='opls_unint_dat.py',dname='opls_unint_dat')
dict=opls_nacl_dat = buildResDict('opls_nacl.in')
write(dict, outfile='opls_nacl_dat.py',dname='opls_nacl_dat')
dict=unict_dat = buildResDict('unict.in')
write(dict, outfile='unict_dat.py',dname='unict_dat')
dict=opls_unict_dat = buildResDict('opls_unict.in')
write(dict, outfile='opls_unict_dat.py',dname='opls_unict_dat')
dict=uni_dat = buildResDict('uni.in')
write(dict, outfile='uni_dat.py',dname='uni_dat')
dict=opls_uni_dat = buildResDict('opls_uni.in')
write(dict, outfile='opls_uni_dat.py',dname='opls_uni_dat')
dict=nh2e_dat = buildResDict('nh2e.in')
write(dict, outfile='nh2e_dat.py',dname='nh2e_dat')
dict=allnt_dat = buildResDict('allnt.in')
write(dict, outfile='allnt_dat.py',dname='allnt_dat')
dict=allct_dat = buildResDict('allct.in')
write(dict, outfile='allct_dat.py',dname='allct_dat')
dict=all_dat = buildResDict('all.in')
write(dict, outfile='all_dat.py',dname='all_dat')
dict=all_nuc94_dat = buildResDict('all_nuc94.in')
write(dict, outfile='all_nuc94_dat.py',dname='all_nuc94_dat')
dict=all_aminoct94_dat = buildResDict('all_aminoct94.in')
#write(dict, outfile='test_all_aminoct94_dat.py',dname='all_aminoct94_dat')
write(dict, outfile='all_aminoct94_dat.py',dname='all_aminoct94_dat')
dict=all_aminont94_dat = buildResDict('all_aminont94.in')
write(dict, outfile='all_aminont94_dat.py',dname='all_aminont94_dat')
dict=all_amino94_dat = buildResDict('all_amino94.in')
write(dict, outfile='all_amino94_dat.py',dname='all_amino94_dat')
