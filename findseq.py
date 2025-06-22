"""
DESCRIPTION:
    Given a sequence/regex to find, select those matching amino acids in the
    protein.

USAGE:
    findseq needle, [haystack[, selName[, het[, matchMode]]]]

PARAMS:
needle (string)
        the sequence of amino acids to match and select
        in the haystack.  This can be a sequence of amino
        acids, or a string-style regular expression.  See
        examples.

haystack (string; PyMOL object or selection; defaults to *)
        name of the PyMOL object/selection in which
        to find the needle.

selName (string; defaults to None)
        This is the name of the selection to return.  If selName
        is left blank (None), then the selection name will be
        foundSeqXYZ where XYZ is some random number; if selName is
        "sele" the usual PyMOL "(sele)" will be used; and, lastly,
        if selName is anything else, that name will be used verbatim.

het (0 or 1; defaults to 0)
        This boolean flag allows (1) or disallows (0) heteroatoms
        from being considered.

matchMode (first/all/chain; defaults to "all")
        Subsequences or motifs might be repeated, this controls how we
        consider multiple matches. Options are:
        - 'first': Return the first match found in each object.
        - 'chain': Return the first match found in each chain.
        - 'all': Return all matches found in all chains.

RETURNS:
    a newly created selection with the atoms you sought.

EXAMPLE:
    # find SPVI in 1h12, foundSeqXYZ as return name
    findseq SPVI, 1h12

    # find FATEW and make it (sele).
    findseq FATEW, 1g01, sele

    # find the regular expression GMS.*QWY in 1a3h
    # and put the return value in (sele).
    fetch 1a3h
    # this ends up finding the sequence, GMSSHGLQWY
    findseq GMS.*QWY, 1a3h, sele
"""

from pymol import cmd
import re
from collections import defaultdict
from itertools import groupby


def findseq(needle, haystack='*', selName=None, het=0, matchMode="all"):
    """
    DESCRIPTION:
        Given a sequence/regex to find, select those matching amino acids in the
        protein.

    USAGE:
        findseq needle, [haystack[, selName[, het[, matchMode]]]]

    PARAMS:
    needle (string)
            the sequence of amino acids to match and select
            in the haystack.  This can be a sequence of amino
            acids, or a string-style regular expression.  See
            examples.

    haystack (string; PyMOL object or selection; defaults to *)
            name of the PyMOL object/selection in which
            to find the needle.

    selName (string; defaults to None)
            This is the name of the selection to return.  If selName
            is left blank (None), then the selection name will be
            foundSeqXYZ where XYZ is some random number; if selName is
            "sele" the usual PyMOL "(sele)" will be used; and, lastly,
            if selName is anything else, that name will be used verbatim.

    het (0 or 1; defaults to 0)
            This boolean flag allows (1) or disallows (0) heteroatoms
            from being considered.

    matchMode (first/all/chain; defaults to "all")
            Subsequences or motifs might be repeated, this controls how we
            consider multiple matches. Options are:
            - 'first': Return the first match found in each object.
            - 'chain': Return the first match found in each chain.
            - 'all': Return all matches found in all chains.

    RETURNS:
        a newly created selection with the atoms you sought.

    EXAMPLE:
        # find SPVI in 1h12, foundSeqXYZ as return name
        findseq SPVI, 1h12

        # find FATEW and make it (sele).
        findseq FATEW, 1g01, sele

        # find the regular expression GMS.*QWY in 1a3h
        # and put the return value in (sele).
        fetch 1a3h
        # this ends up finding the sequence, GMSSHGLQWY
        findseq GMS.*QWY, 1a3h, sele
    """
    selName = selName or cmd.get_unused_name("foundSeq")

    # input checking
    if not checkParams(needle, haystack, selName, het, matchMode):
        print("There was an error with a parameter. Please see the above error message for how to fix it.")
        return

    # remove hetero atoms (waters/ligands/etc) from consideration?
    haystack_sel = cmd.get_unused_name()
    if bool(int(het)):
        cmd.select(haystack_sel, f"({haystack})")
    else:
        cmd.select(haystack_sel, f"({haystack}) and not het")

    # get the AAs in the haystack
    data = defaultdict(lambda: defaultdict(list))
    for obj in cmd.get_object_list():
        cmd.iterate(f"%{obj} and {haystack_sel}", lambda atom: data[atom.model][(atom.segi, atom.chain)].append((atom.resi, atom.oneletter)))

    reNeedle = re.compile(needle.upper())
    matches = []
    for obj, val in data.items():
        found = False
        for (segi, chain), seq_data in val.items():
            # discard repeated atoms per residue
            seq = "".join(resn for (resi, resn), _ in groupby(seq_data))
            ids = [resi for (resi, resn), _ in groupby(seq_data)]

            for m in reNeedle.finditer(seq):
                found = True
                start, stop = m.span()
                resi = "+".join(ids[start:stop])

                sel = f"/{obj}/{segi}/{chain}/{resi}"
                matches.append(sel)

                if matchMode == "chain" or matchMode == "first":
                    break

            if matchMode == "first" and found:
                break

    # select all residues that match
    cmd.delete(haystack_sel)
    if not matches:
        print("Sequence was not found")
        return

    cmd.select(selName, " or ".join(matches))
    print(f'findseq: selection "{selName}" defined with {cmd.count_atoms(selName)} atoms.')
    return selName


cmd.extend("findseq", findseq)
cmd.findseq = findseq
cmd.auto_arg[1]['findseq'] = [ cmd.object_sc, 'object', '']
cmd.auto_arg[2]['findseq'] = [lambda: cmd.Shortcut(['het=1','matchMode=chain', "matchMode=first"]), 'params', '']
cmd.auto_arg[3]['findseq'] = [lambda: cmd.Shortcut(['het=1','matchMode=chain', "matchMode=first"]), 'params', '']


def checkParams(needle, haystack, selName, het, matchMode):
    """
    This is just a helper function for checking the user input
    """
    # check Needle
    if not isinstance(needle, str) or len(needle) == 0:
        print("Error: Please provide a string 'needle' to search for.")
        return False

    # check Haystack
    if not isinstance(haystack, str) or len(haystack) == 0:
        print("Error: Please provide valid PyMOL object or selection name in which to search.")
        return False

    # check het
    try:
        het = bool(int(het))
    except ValueError:
        print("Error: The 'het' parameter was not 0 or 1.")
        return False

    # check matchMode
    if matchMode not in ("first", "chain", "all"):
        print("Error: 'matchMode' parameter must be one of 'first', 'chain' or 'all'")
        return False

    # check selName
    if not isinstance(selName, str):
        print("Error: selName was not a string.")
        return False
    return True


def test_findseq():
    cmd.reinitialize()
    cmd.fab('ACDEFG', 'm1')
    cmd.fab('HIKLMN', 'm2')
    cmd.do('findseq DEF, *, foundSeq')
    assert cmd.count_atoms("foundSeq and m1") == 47
    assert cmd.count_atoms("foundSeq and m2") == 0
