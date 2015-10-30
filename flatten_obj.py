from pymol import cmd, stored
import re
try:
    from collections import OrderedDict
    _orderedDict = True
except ImportError:
    _orderedDict = False

# PyMOL 1.7.4 introduces support for multi-letter chains, so we can afford to
# use a smaller alphabet. In earlier versions, use lower-case letters if needed
# (requires running `set ignore_case, 0`)
_long_chains = cmd.get_version()[1] >= 1.74
_default_base = 36 if _long_chains else 62

class OutOfChainsError(Exception):
    def __init__(self,msg):
        self.msg=msg
    def __str__(self):
        return str(self.msg)

class ChainSet(object):
    """
    Base class for various methods to rename chains

    Contains _chains, which maps from the renamed chain to a tuple with the
    original (object,state,chain). All dict-like accessors work on ChainSets,
    e.g.
        chain_set["A"] -> ("obj",1,"A")

    """
    def __init__(self):
        # Use an OrderedDict in Python >= 1.7 for better printing
        if _orderedDict:
            self._chains = OrderedDict()
        else:
            self._chains = dict()

    def map_chain(self, obj, state, origChain ):
        """
        map_chain(string obj,int state, string chain]]) -> string

        Maps a chain letter to a unique chainID. Results are unique within each
        instance, and can be used as keys on this chain set.
        """
        raise NotImplementedError("Base class")

    # delegate most methods to _chains
    def __getattr__(self,at):
        if at in "pop popitem update setdefault".split():
            raise AttributeError("type object '%s' has no attribute '%s'"%(type(self),at))
        return getattr(self._chains,at)
    def __cmp__(self,other):    return self._chains.__cmp__(other)
    def __eq__(self,other):     return self._chains.__eq__(other)
    def __ge__(self,other):     return self._chains.__ge__(other)
    def __gt__(self,other):     return self._chains.__gt__(other)
    def __le__(self,other):     return self._chains.__le__(other)
    def __lt__(self,other):     return self._chains.__lt__(other)
    def __ne__(self,other):     return self._chains.__ne__(other)
    def __len__(self):          return self._chains.__len__()
    def __contains__(self,key): return self._chains.__contains__(key)
    def __getitem__(self,key):  return self._chains.__getitem__(key)
    def __iter__(self):         return self._chains.__iter__()
    def __str__(self):          return str(self._chains)

    @staticmethod
    def _int_to_chain(i,base=_default_base):
        """
        _int_to_chain(int,int) -> str

        Converts a positive integer to a chain ID. Chain IDs include uppercase
        characters, numbers, and optionally lowercase letters.

        i = a positive integer to convert
        base = the alphabet size to include. Typically 36 or 62.
        """
        if i < 0:
            raise ValueError("positive integers only")
        if base < 0 or 62 < base:
            raise ValueError("Invalid base")

        quot = int(i)/base
        rem = i%base
        if rem < 26:
            letter = chr( ord("A") + rem)
        elif rem < 36:
            letter = str( rem-26)
        else:
            letter = chr( ord("a") + rem - 36)
        if quot == 0:
            return letter
        else:
            return ChainSet._int_to_chain(quot-1,base) + letter


class DefaultChainSet(ChainSet):
    """
    Avoids relettering chains if possible. If a chain has been used, uses the
    next available chain letter. Note that this can potentially lead to
    cascading renames, e.g. if chains are sorted alphabetically rather than by
    object.

    Used for rename = 0.
    """
    def __init__(self):
        super(DefaultChainSet,self).__init__()
        self._next_chain = 0
    def map_chain(self, obj, state, origChain ):
        # Keep _next_chain up-to-date
        while ChainSet._int_to_chain(self._next_chain) in self:
            self._next_chain += 1
        # Map this chain
        if origChain in self:
            # Rename
            next_chain = ChainSet._int_to_chain(self._next_chain)
            self._next_chain += 1
        else:
            next_chain = origChain
        self._chains[next_chain] = (obj,state,origChain)
        return next_chain

class SequentialChainSet(ChainSet):
    """
    Renumbers all chains starting at A, continuing through the capital letters
    and numbers, and then adding additional letters through 9999 (the last
    valid chain for mmCIF) and beyond.

    Used for rename=1
    """
    def __init__(self):
        super(SequentialChainSet,self).__init__()
        self._next_chain = 0

    def map_chain(self, obj, state, origChain ):
        next_chain = ChainSet._int_to_chain(self._next_chain)
        self._chains[next_chain] = (obj,state,origChain)
        self._next_chain += 1
        return next_chain

class LongChainSet(ChainSet):
    """
    Uses long strings for the chain names. Chains are renamed like
    "%s_%s_%04d"%(original_chainid,objectname,state).

    Used for rename=2
    """
    def map_chain(self, obj, state, origChain ):
        ch = "%s_%s_%04d"%(origChain,obj,state)
        if self.has_key(ch):
            raise ValueError("Duplicate chain %s"%(ch))
        self._chains[ch] = (obj,state,origChain)
        return ch




def flatten_obj(name="",selection="",state=0,rename=0,quiet=1,chain_map=""):
    """
DESCRIPTION

    "flatten_obj" combines multiple objects or states into a single object,
    renaming chains where required

USAGE

    flatten_obj name, selection[, state[, rename[, quiet[, chain_map]]]]

ARGUMENTS

    name = a unique name for the flattened object {default: flat}

    selection = the set of objects to include in the flattening. The selection
        will be expanded to include all atoms of objects. {default: all}

    state = the source state to select. Use 0 or -1 to flatten all states {default: 0}

    rename = The scheme to use for renaming chains: {default: 0}
        (0) preserve chains IDs where possible, rename other chains
            alphabetically
        (1) rename all chains alphabetically
        (2) rename chains using the original chain letter, object name, and state

    quiet = If set to 0, print some additional information about progress and
        chain renaming {default: 1}

    chain_map = An attribute name for the 'stored' scratch object. If
        specified, `stored.<chain_map>` will be populated with a dictionary
        mapping the new chain names to a tuple giving the originated object,
        state, and chainID. {default: ""}

NOTES

    Like the select command, if name is omitted then the default object name
    ("flat") is used as the name argument.

    Chain renaming is tricky. PDB files originally limited chains to single
    letter identifiers containing [A-Za-z0-9]. When this was found to be
    limiting, multi-letter chains (ideally < 4 chars) were allowed. This is
    supported as of PyMOL 1.7. Earlier versions do not accept rename=2, and
    will raise an exception when flattening a structure with more than 62
    chains.

EXAMPLES

    flatten_obj flat, nmrObj
    flatten_obj ( obj1 or obj2 )

PYMOL API

    cmd.flatten_obj(string name, string selection)

SEE ALSO

    split_states

    """

    # arguments

    # Single argument; treat as selection
    if name and not selection:
        selection = name
        name = ""
    # default name and selection
    if not name:
        name = "flat"
    if not selection:
        selection = "(all)"

    state = int(state)
    rename = int(rename)
    quiet = int(quiet)

    # Wrap in extra parantheses for get_object_list
    selection = "( %s )" % selection

    if rename == 0:
        chainSet = DefaultChainSet()
    elif rename == 1:
        chainSet = SequentialChainSet()
    elif rename == 2:
        chainSet = LongChainSet()
    else:
        raise ValueError("Unrecognized rename option (Valid: 0,1,2)")

    metaprefix = "temp" #TODO unique prefix

    # store original value of retain_order, which causes weird interleaving of
    # structures if enabled.
    retain_order = cmd.get("retain_order")
    try:
        cmd.set("retain_order",0)

        # create new object for each state
        for obj in cmd.get_object_list(selection):

            if state <= 0:
                # all states
                prefix = "%s_%s_"%(metaprefix,obj)
                cmd.split_states(obj,prefix=prefix)
            else:
                prefix = "%s_%s_%04d"%(metaprefix,obj,state)
                cmd.create(prefix, obj, state, 1)

        # renumber all states
        statere = re.compile("^%s_(.*)_(\d+)$" % metaprefix) # matches split object names

        warn_lowercase = False

        # Iterate over all objects with metaprefix
        try:
            for obj in cmd.get_object_list("(%s_*)"%(metaprefix) ):
                m = statere.match(obj)
                if m is None:
                    print("Failed to match object %s" %obj)
                    continue
                origobj = m.group(1)
                statenum = int(m.group(2))

                chains = set([atom.chain for atom in cmd.get_model(obj).atom])

                rev_chain_map = {} #old -> new, for this obj only
                for chain in sorted(chains):
                    new_chain = chainSet.map_chain(origobj,statenum,chain)
                    rev_chain_map[chain] = new_chain
                    if not quiet:
                        print("  %s state %d chain %s -> %s"%(origobj,statenum,chain, new_chain) )
                    if not _long_chains:
                        if len(new_chain) > 1:
                            raise OutOfChainsError("No additional chains available (max 62).")

                try:
                    stored_name = stored.get_unused_name()
                except AttributeError:
                    # backup solution for PyMOL < 1.6
                    stored_name = "_flatten_obj_chainMapping"

                setattr(stored,stored_name,rev_chain_map)
                cmd.alter(obj,"chain = stored.%s[chain]"%stored_name)
                delattr(stored,stored_name)

            print("Creating object from %s_*"%metaprefix)
            # Recombine into a single object
            cmd.create(name,"%s_*"%metaprefix)

            # Set chain_map
            if chain_map:
                setattr(stored,chain_map,chainSet)

            # Warn if lowercase chains were generated 
            if cmd.get("ignore_case") == "on" and any([c.upper() != c for c in chainSet.keys()]):
                print("Warning: using lower-case chain IDs. Consider running the "
                        "following command:\n  set ignore_case, 0" )

        finally:
            # Clean up
            print("Cleaning up intermediates")
            cmd.delete("%s_*"%metaprefix)
    finally:
        # restore original parameters
        print("Resetting variables")
        cmd.set("retain_order",retain_order)


cmd.extend('flatten_obj', flatten_obj)

# tab-completion of arguments
#cmd.auto_arg[3]['flatten_obj'] = [ cmd.object_sc, 'object', '']
