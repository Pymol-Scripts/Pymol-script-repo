# very rough bounds!
typemap = { (-180,0, 90, 180) : "BR", (-150, -30, -60, 60) : "AR", (0, 180, 90, 180) : "BL", (30, 150, -60, 60) : "AL" }
 
def determinetype(phipsi):
    phi, psi = phipsi
    for bound in typemap:
        if bound[0] < phi < bound[1] and bound[2] < psi < bound[3]:
            return typemap[bound]
    return "?"
 
def my_phi_psi(selection):
    r = cmd.get_phipsi(selection)
 
    if r is not None:
        keys = r.keys()
        keys.sort()
 
        cmd.feedback('push')
        cmd.feedback('disable','executive','actions')
        for key in keys:
            phipsiType = determinetype(r[key])
            argtuple = r[key] + (phipsiType,)
            cmd.iterate("(%s`%d)" % key, "print ' %-5s " + ("( %4.0f, %4.0f ) %s" % argtuple) + "'%(resn+'-'+resi+' '+chain+':')")
        cmd.feedback('pop')
 
def motif(residueRange, chain=None):
    """
    Use like "motif 10-12" to show a backbone segment from residues 10 to 12.
    """
    name = residueRange
    if chain is None:
        selection = "name ca+n+c+o+h and not hetatm and resi %s" % residueRange
    else:
        selection = "name ca+n+c+o+h and not hetatm and resi %s and chain %s" % (residueRange, chain)
    cmd.select(name, selection)
    cmd.show("sticks", name)
    cmd.zoom(name)
    cmd.disable(name)
    cmd.label("name ca and %s" % selection, "resi")
    my_phi_psi(name)
 
cmd.extend("motif", motif)
cmd.extend("my_phi_psi", my_phi_psi)

