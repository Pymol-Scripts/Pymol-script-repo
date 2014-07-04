from pymol import cmd


def removealt(obj="(all)", keep="A"):
    """
    removeAlt -- remove all alternate location-atoms not of altloc "keep" from object.

    input:
            obj -- the object(s) to remove the atoms frmo
            keep -- which type of alt loc to keep

    output: none -- removes atoms

    examples:
            removeAlt # remove all altLocations that aren't altloc A
            removeAlt pdbID, C  # remove all but C altlocations from pdbID
    """
    # select & remove all non A altlocs
    remStr = "%s and not (alt ''+%s)" % (obj, keep)
    cmd.remove(remStr)
    # reset the PDB information
    cmd.alter(obj, "alt=''")

cmd.extend("removealt", removealt)
