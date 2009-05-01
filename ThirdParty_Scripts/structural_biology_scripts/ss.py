from pymol import cmd
from pymol import stored
 
def ss(selection):
 
    class SSE(object):
 
        def __init__(self, start, typ):
            self.start, self.typ = start, typ
            self.end = -1
 
        def __repr__(self):
            return "%s-%s %s" % (self.start, self.end, self.typ)
 
    stored.pairs = []
    cmd.iterate("%s and n. ca" % selection, "stored.pairs.append((resi, ss))")
    num, currentType = stored.pairs[0]
 
    sses = [SSE(num, currentType)]
    currentSSE = sses[0]
    for resi, ssType in stored.pairs:
        if ssType == currentType:
            currentSSE.end = resi
        else:
            sses.append(SSE(resi, ssType))
            currentSSE = sses[-1]
            currentType = ssType
 
    for sse in sses:
        print sse
 
cmd.extend("ss", ss)

