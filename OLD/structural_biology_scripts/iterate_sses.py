def iterate_sses(selection, action):
 
    class SSE(object):
 
        def __init__(self, start, typ, sseNumber):
            self.start, self.typ = start, typ
            self.end = -1
            self.sseNumber = sseNumber
 
        def __repr__(self):
            return "%s-%s %s" % (self.start, self.end, self.typ)
 
    stored.pairs = []
    cmd.iterate(selection, "stored.pairs.append((resi, ss))")
    num, currentType = stored.pairs[0]
 
    sseNumber = 1
    sses = [SSE(num, currentType, sseNumber)]
    currentSSE = sses[0]
    for resi, ssType in stored.pairs:
        if ssType == currentType:
            currentSSE.end = resi
        else:
            sseNumber += 1
            sses.append(SSE(resi, ssType, sseNumber))
            currentSSE = sses[-1]
            currentType = ssType
 
    for sse in sses:
        action(sse)
 
cmd.extend("iterate_sses", iterate_sses)
