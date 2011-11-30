from pymol import cmd
from itertools import izip
import math
 
def rgyrate(selection='(all)', quiet=1):
    '''
DESCRIPTION
 
    Radius of gyration
 
USAGE
 
    rgyrate [ selection ]
    '''
    quiet = int(quiet)
    model = cmd.get_model(selection).atom
    x = [i.coord for i in model]
    mass = [i.get_mass() for i in model]
    xm = [(m*i,m*j,m*k) for (i,j,k),m in izip(x,mass)]
    tmass = sum(mass)
    rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in izip(x,xm))
    mm = sum((sum(i)/tmass)**2 for i in izip(*xm))
    rg = math.sqrt(rr/tmass - mm)
    if not quiet:
        print "Radius of gyration: %.2f" % (rg)
    return rg
 
cmd.extend("rgyrate", rgyrate)

