from __future__ import print_function
from pymol import cmd


def findSurfaceAtoms(selection="all",cutoff=2.5):
    """
    Adapted from Jason Vertrees https://pymolwiki.org/index.php/FindSurfaceResidues
DESCRIPTION

    Finds those atoms on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.

USAGE

    findSurfaceAtoms [ selection, [ cutoff ]]

SEE ALSO

    findSurfaceResidues
    """
    cutoff = float(cutoff)

    tmpObj = cmd.get_unused_name("_tmp")
    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)

    cmd.set("dot_solvent", 1, tmpObj)
    cmd.get_area(selection=tmpObj, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    cmd.remove(tmpObj + " and b < " + str(cutoff))

    selName = cmd.get_unused_name("exposed_atm_")

    cmd.select(selName, "(" + selection + ") in " + tmpObj)

    cmd.delete(tmpObj)

    return selName


def findSurfaceCharge(pH=7.0, folded=True, selection="all", cutoff=2.5):
    """
DESCRIPTION

    Calculates a surface charge at entered pH. Also allows for the charge of an unfolded protein to be calculated.

USAGE

    findSurfaceCharge [pH, [folded, [selection ,[cutoff]]]]

ARGUMENTS

    pH = The pH value to estimate a surface charge at

    folded = Whether the protein is folded (True) or denatured (False)

    selection = string: object or selection in which to find exposed
    residues {default: all}

    cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

RETURNS

    A printout of the estimated surface charge at a given pH

    """
    cutoff = float(cutoff)

    selName = findSurfaceAtoms(selection, cutoff)

    exposed = set()
    cmd.iterate(selName, "exposed.add((resv))", space=locals())
    cmd.delete(selName)

    selNameRes = cmd.get_unused_name("exposed_res_")

    exposed=sorted(exposed) #list of exposed residues

    seq=cmd.get_fastastr('all')
    seqbegin = seq.find('\n')
    newSeq = seq[seqbegin::].replace('\n', ' ').replace('\r', '').replace(' ','').replace("?","")

    #adjusts for beginning position
    first = set()
    allRes = findSurfaceAtoms(selection, 0)
    cmd.iterate(allRes, "first.add((resv))", space=locals())
    cmd.delete(allRes)

    selNameRes = cmd.get_unused_name("exposed_res_")

    first=sorted(first)[0] #firstRes
    #gets all charged amino acids on the surface
    reslist= []
    exposedAtms=""
    K=0
    R=0
    D=0
    H=0
    E=0

    if folded:
        offset=(1+first)
        for r in exposed:
            amino=newSeq[r-offset]
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
            exposedAtms+=amino
            chargedAA=amino
            reslist.append(chargedAA)
    else:
        for r in newSeq:
            amino=r
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
            exposedAtms+=amino
            chargedAA=amino
            reslist.append(chargedAA)


    kCharge= 1 / (1 + 10 ** (pH - 10.54))
    rCharge= 1 / (1 + 10 ** (pH - 12.48))
    dCharge= -(1 / (1 + 10 ** (4.07 - pH)))
    eCharge= -(1 / (1 + 10 ** (3.90 - pH)))
    hCharge= 1 / (1 + 10 ** (pH - 6.04))

    charge=kCharge*K+rCharge*R+hCharge*H+dCharge*D+eCharge*E
    charge=round(charge,2)
    if charge>0:
        chargetx="+"+str(charge)
    else:
        chargetx=str(charge)

    if folded:
        print ("Exposed charged residues: " +str(exposedAtms))
        print ("The expected surface charge of this protein at pH " + str(pH) +" is: " +chargetx)

    else:
        print ("Charged residues: "+str(exposedAtms))
        print ("The expected charge of this denatured protein at pH " +str(pH) +" is: " +chargetx)


cmd.extend("findSurfaceAtoms", findSurfaceAtoms)
cmd.extend("findSurfaceCharge", findSurfaceCharge)
