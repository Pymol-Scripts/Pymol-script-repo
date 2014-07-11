import pymol
from pymol import cmd


def toGroup(groupName, sel, prefix="", delOrig=True):
    """
    DESCRIPTION
    toGroup will take a multistate object and extract it
    to a group with N objects all in state #1.  It essentially
    performs the following:

    split_states myObj, prefix=somePrefix
    group newGroup, somePrefix*
    delete myObj

    PARAMETERS:

    groupName (string)
        The name of the group to create

    sel (string)
        The name of the selection/object from which
        to make the group

    prefix (string)
        The prefix of the names of each of split states.
        For example, if your prefix is ''obj'' and is in
        states 1 through 100 then the states will be labeled
        obj1, obj2, obj3, ..., obj100.

    delOrig (string/boolean)
        If true then delete the original selection, otherwise not.

    RETURN

    Nothing, it makes a new group.

    """
    if prefix == "":
        prefix = sel + "_grouped"

    cmd.split_states(sel, prefix=prefix)
    cmd.group(groupName, prefix + "*")

    if delOrig:
        cmd.delete(sel)

cmd.extend("toGroup", toGroup)
