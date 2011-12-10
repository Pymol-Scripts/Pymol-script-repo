## def uniq(l, func=None):
##     """Return a new list with duplicate items removed."""

##     l2 = l[:]	# make a copy
##     d = {}
##     def add_to_dict(value,d=d):
## 	d[`value`] = value
##     map(add_to_dict,l2)
##     l3 = d.values()
##     if len(l2)==len(l3): return(l2)
##     if func: l3.sort(func)
##     else: l3.sort()
##     return l3

def uniq(alist):    # Fastest order preserving
    set = {}
    return [set.setdefault(e,e) for e in alist if e not in set]
 
def uniq3(alist):    # Fastest without order preserving
    set = {}
    map(set.__setitem__, alist, [])
    return set.keys()

"""
from mglutil.util.uniq import uniq, uniq2, uniq3
import time
a=range(100)
b=range(10)
c=a+b

t1=time.time()
for i in range(5000):  x=uniq(c)
print time.time()-t1

t1=time.time()
for i in range(5000):  x=uniq2(c)
print time.time()-t1

t1=time.time()
for i in range(5000):  x=uniq3(c)
print time.time()-t1


>>>
0.865363121033
0.463307857513
0.260641098022
"""
