'''
http://pymolwiki.org/index.php/FindSurfaceResiduesListCharged
'''

from __future__ import print_function
from pymol import cmd


def findSurfaceAtoms(selection="all",cutoff=2.5, quiet=1):
    """
DESCRIPTION

    Finds those atoms on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.

USAGE

    findSurfaceAtoms [ selection, [ cutoff ]]

SEE ALSO

    findSurfaceResidues
    """
    cutoff, quiet = float(cutoff), int(quiet)

    tmpObj = cmd.get_unused_name("_tmp")
    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)

    cmd.set("dot_solvent", 1, tmpObj)
    cmd.get_area(selection=tmpObj, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    cmd.remove(tmpObj + " and b < " + str(cutoff))

    selName = cmd.get_unused_name("exposed_atm_")
    cmd.select(selName, "(" + selection + ") in " + tmpObj)

    cmd.delete(tmpObj)

    if not quiet:
        print("Exposed atoms are selected in: " + selName)

    return selName


def findSurfaceResiduesListCharged(pH=7.0, selection="all", cutoff=2.5, doShow=0, quiet=1):
    """
DESCRIPTION

    Identifies and lists all charged surface residues. Also calculates a
    surface charge at entered pH.

USAGE

    findSurfaceResiduesListCharged [pH, [selection ,[cutoff ,[ doShow , [ quiet]]]]]]

ARGUMENTS

    pH = The pH value to estimate a surface charge at

    selection = string: object or selection in which to find exposed
    residues {default: all}

    cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

RETURNS

    (list: (chain, resv ) )
        A Python list of residue numbers corresponding
        to those residues w/more exposure than the cutoff.

    A printout of all charged amino acids and the estimated surface charge of a protein.

    """
    cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)

    selName = findSurfaceAtoms(selection, cutoff, quiet)

    exposed = set()
    cmd.iterate(selName, "exposed.add((chain,resv))", space=locals())

    selNameRes = cmd.get_unused_name("exposed_res_")
    cmd.select(selNameRes, "byres " + selName)

    if not quiet:
        print("Exposed residues are selected in: " + selNameRes)

    if doShow:
        cmd.show_as("spheres", "(" + selection + ") and polymer")
        cmd.color("white", selection)
        cmd.color("yellow", selNameRes)
        cmd.color("red", selName)
    exposed=sorted(exposed) #list of exposed residues
    seq=cmd.get_fastastr('all')
    newseq = seq.replace('\n', ' ').replace('\r', '').replace(' ','').replace("?","")
    nDomains=newseq.count(">")

    def splitSeq(seq,d): #splits residues into their domains
        domain=-1
        count=0
        splits=[""]*d
        for r in seq:
            if r==">":
                domain+=1
                count=7
            if count>0:
                count-=1
                continue
            if r=="\n":
                continue
            splits[domain]+=r
        return splits

    splits=splitSeq(newseq,nDomains) #list of all residues separated by domain

    #gets all charged amino acids on the surface
    reslist=[[] for x in range(nDomains)]
    exposedAms=""
    K=0
    R=0
    D=0
    H=0
    E=0

    for r in exposed:
        tempDomain=ord(r[0])-65
        if r[1]-1 >= len(splits[tempDomain]):
            continue
        amino=splits[tempDomain][r[1]-1]
        if amino not in "KRDHE":
            continue
        elif amino=='K':
            K+=1
        elif amino=='R':
            R+=1
        elif amino=='D':
            D+=1
        elif amino=="H":
            H+=1
        elif amino=='E':
            E+=1
        exposedAms+=amino
        amPos=amino+str(r[1])
        reslist[tempDomain].append(amPos)

    n=1
    for dom in reslist:
        print("Domain " + str(n)+" charged residues:")
        print(dom)
        n+=1

    print("Charged residue list")
    print(exposedAms)
    print("Number of exposed charged residues")
    print(len(exposedAms))

    if pH > 10.54:
        kCharge=10 ** -(pH-10.54)
    else:
        kCharge=1-(10** (pH-10.54))

    if pH > 12.48:
        rCharge=10 ** -(pH-12.48)
    else:
        rCharge=1-(10** (pH-12.48))

    if pH > 4.07:
        dCharge=-(1-(10**-(pH-4.07)))
    else:
        dCharge=-(10**(pH-4.07))

    if pH > 3.90:
        eCharge=-(1-(10**-(pH-3.90)))
    else:
        eCharge=-(10**(pH-3.90))

    if pH >6.04:
        hCharge=10 ** -(pH-6.04)
    else:
        hCharge=1-(10** (pH-6.04))


    charge=kCharge*K+rCharge*R+hCharge*H+dCharge*D+eCharge*E
    charge=round(charge,2)
    if charge>0:
        chargetx="+"+str(charge)
    else:
        chargetx=str(charge)


    print ("The expected surface charge of this protein at pH " + str(pH) +" is: " +chargetx)
    return sorted(exposed)

cmd.extend("findSurfaceAtoms", findSurfaceAtoms)
cmd.extend("findSurfaceResiduesListCharged", findSurfaceResiduesListCharged)
