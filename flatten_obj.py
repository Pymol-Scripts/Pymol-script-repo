from pymol import cmd, stored
import re

def flatten_obj(*args,**kwargs):
    """
    flatten_obj(name,selection)

DESCRIPTION

    "flatten_obj" combines multiple objects or states into a single object,
    renaming chains where required

USAGE

    flatten_obj name, selection

ARGUMENTS

    name = a unique name for the flattened object

    selection = the set of objects to include in the flattening. The selection
        will be expanded to include all atoms of objects.

NOTES

    Like the select command, if a selection-expression with explicit
    surrounding parethenses is provided as the first argument, then the default
    object name ("flattened") is used as the name argument.

EXAMPLES

    flatten_obj flattened, nmrObj
    flatten_obj ( obj1 or obj2 )

PYMOL API

    cmd.flatten_obj(string name, string selection)

SEE ALSO

    split_states

    """

    # arguments
    name = kwargs.get("name") or "flattened"
    selection = kwargs.get("selection") or "(all)"

    if len(args) == 1:
        # check for paranthetic string
        first = args[0].strip()
        if first[0] == '(' and first[-1] == ')':
            selection = args[0]
        else:
            name = args[0]
    elif len(args) > 1:
        name = args[0]
    if len(args) >= 2:
        name = args[0]
        selection = args[1]

    #print name,selection

    # Wrap in extra parantheses for get_object_list
    selection = "( %s )" % selection


    nextChain = "A"

    metaprefix = "temp" #unique prefix

    # create new object for each state
    for object in cmd.get_object_list(selection):
        prefix = "%s_%s_"%(metaprefix,object)

        cmd.split_states(object,prefix=prefix)

    # renumber all states
    usedstates = set()
    statere = re.compile("^%s_(.*)_(\d+)$" % metaprefix) # matches split object names
    print("objects: "+repr(cmd.get_object_list("%s_*"%(metaprefix) )))

    # Iterate over all objects with metaprefix
    # Globs work sporadically with get_object_list, so iterate all
    for object in cmd.get_object_list("(%s_*)"%(metaprefix) ):
        m = statere.match(object)
        if m is None:
            print("Conflict with object name %s" %object)
            continue
        origobj, statenum = m.groups()

        chains = set([atom.chain for atom in cmd.get_model(object).atom])
        print("Chains for %s: %s" %(object,repr(chains)))
        chainMapping = {}
        for chain in chains:

            if chain in usedstates:
                # Find next unused chain
                while nextChain in usedstates:
                    nextChain = chr(ord(nextChain)+1)

                if nextChain > "Z":
                    raise RuntimeError("No additional chains available!")

                #print( "Renaming chain %s to %s" % (chain, nextChain))

                chainMapping[chain] = nextChain
                usedstates.add(nextChain)
            else:
                #print( "Chain %s free" % chain)
                chainMapping[chain] = chain
                usedstates.add(chain)
        stored._flatten_obj_chainMapping = chainMapping

        for c1,c2 in chainMapping.items():
            if c1 != c2:
                print("Renaming %s state %s chain %s to chain %s" % (origobj, statenum, c1, c2) )

        cmd.alter(object,"chain = stored._flatten_obj_chainMapping[chain]")

        #print( "usedstates="+repr(usedstates))

    # Recombine into a single object
    cmd.create(name,"%s_*"%metaprefix)

    # Clean up
    cmd.delete("(%s_*)"%metaprefix)
#   nextchain = "A"
#     for state in xrange(1):
#         atoms = cmd.get_model(selection)
#         for atom in atoms.atom:
#             print atom.chain
#             atom.chain = "Z"
#             print atom.chain


cmd.extend('flatten_obj', flatten_obj)

# tab-completion of arguments
#cmd.auto_arg[3]['flatten_obj'] = [ cmd.object_sc, 'object', '']
