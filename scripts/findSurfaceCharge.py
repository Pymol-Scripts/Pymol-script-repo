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


def _findSurfaceChargeImpl(selection, pH, folded, cutoff):

    def get_exposed_residues(selection,cutoff):
        cutoff = float(cutoff)


        selName = findSurfaceAtoms(selection, cutoff)

        tempExposed = set()
        cmd.iterate(selName, "tempExposed.add((model,segi,chain,resv,resi,oneletter))", space=locals())
        cmd.delete(selName)

        tempExposed=sorted(tempExposed) #list of exposed residues
        exposed=[]
        for res in tempExposed:
            exposed.append(res[-1] + res[-2])

        return exposed

    if folded:
        exposed = get_exposed_residues(selection,cutoff)
    else:
        exposed = get_exposed_residues(selection,0)

    pH=float(pH)

    #gets all charged amino acids on the surface
    exposedAtms=""
    K=0
    R=0
    D=0
    H=0
    E=0

    for r in exposed:
        amino=r[0]
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


    kCharge= 1 / (1 + 10 ** (pH - 10.54))
    rCharge= 1 / (1 + 10 ** (pH - 12.48))
    dCharge= -(1 / (1 + 10 ** (4.07 - pH)))
    eCharge= -(1 / (1 + 10 ** (3.90 - pH)))
    hCharge= 1 / (1 + 10 ** (pH - 6.04))

    charge=kCharge*K+rCharge*R+hCharge*H+dCharge*D+eCharge*E
    chargetx = "%+.2f" % (charge)

    if folded:
        print ("Exposed charged residues: " +str(exposedAtms))
        print ("The expected surface charge of " + selection +" at pH " + str(pH) +" is: " +chargetx)

    else:
        print ("Charged residues: "+str(exposedAtms))
        print ("The expected charge of denatured " + selection +" at pH " +str(pH) +" is: " +chargetx)
    return (selection, chargetx)


def findSurfaceCharge(selection="", pH=7.0, folded=True, cutoff=2.5):
    """
DESCRIPTION

    Calculates a surface charge at entered pH. Also allows for the charge of an unfolded protein to be calculated.

USAGE

    findSurfaceCharge [pH, [folded, [selection ,[cutoff]]]]

ARGUMENTS

    pH = The pH value to estimate a surface charge at

    folded = Whether the protein is folded (True) or denatured (False)

    selection = string: object or selection in which to find exposed
        residues {default: empty string - all objects}

    cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

RETURNS

    A printout of the estimated surface charge at a given pH

    """
    if not selection:
        for obj in cmd.get_names():
            _findSurfaceChargeImpl(obj, pH, folded, cutoff)
    else:
        _findSurfaceChargeImpl(selection, pH, folded, cutoff)


cmd.extend("findSurfaceAtoms", findSurfaceAtoms)
cmd.extend("findSurfaceCharge", findSurfaceCharge)
