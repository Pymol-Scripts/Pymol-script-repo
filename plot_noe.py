'''
http://pymolwiki.org/index.php/plot_noe

(c) 2011-09-23 by Justin L Lorieau
(c) 2012-10-31 by Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException
import re

def plot_noe(filename, selection='', line_color='gray20', line_width='1.0', single=0, quiet=1, aria=0, per_atom=0, per_residue=1):
    """
DESCRIPTION

    A function for plotting XPLOR NOE restraints on a structure

ARGUMENTS

    filename = string: The filename of the NOE retraint file in XPLOR NIH format.

    selection = string: atom selection {default: all}

    line_color = string: The color for the NOE lines. {default: black}

    line_width = float: The thickness of the NOE lines. {default: 1.0}

    aria = integer: Name NOEs after Aria IDs.

    per_atom: Group NOEs on atom basis

    per_residue: Group NOEs on residue basis (default)

NOE Restraint Format Example

    assign (residue 5 and name HB#) (residue 21 and name HA) 3.0 0.7 0.7

EXAMPLE

    PyMOL> plot_noe noe_short.tbl
    """
    from pymol.parsing import split

    single, quiet, aria, per_residue, per_atom = int(single), int(quiet), int(aria), int(per_residue), int(per_atom)

    count = 0
    parent_id = None

    rgx_restraint = re.compile(
        r"""
                            \s*(?P<parent>\w*)\s+\(
                            (segid\s+\"\s*(?P<src_segid>[\w\d]+)\"\s+and\s+)?
                            resid\s+(?P<src_resid>\d+)\s+and\s+
                            name\s+(?P<src_atom>[\w\d]+)
                            \)\s+\(
                            (segid\s+\"\s+(?P<dst_segid>[\w\d]+)\"\s+and\s+)?
                            resid\s*(?P<dst_resid>\d+)\s+and\s+
                            name\s*(?P<dst_atom>[\w\d]+)
                            \)
                            (\s+(?P<dist>[\d\.]+))?
                            (.*id=(?P<ID>\d+))?
                            .*
                        """, re.X)


    for line in open(filename):
        restraint = rgx_restraint.match(line)
        count += 1

        if restraint:
            restraint = restraint.groupdict()

            if restraint["src_segid"] == None:
                restraint["src_segid"] = "A"

            if restraint["dst_segid"] == None:
                restraint["dst_segid"] = "A"

            if restraint["parent"][:4].upper() == "ASSI":
                parent = restraint

            sele1 = 'segid %s and resi %s and name %s' % (restraint["src_segid"], restraint["src_resid"], restraint["src_atom"])
            sele2 = 'segid %s and resi %s and name %s' % (restraint["dst_segid"], restraint["dst_resid"], restraint["dst_atom"])

            if single:
                label = "NOE"

            elif aria and restraint["ID"]:
                label = "NOE_" + restraint["ID"]
            elif aria and parent and parent["ID"]:
                label = "NOE_" + parent["ID"]

            elif per_atom:
                label = "NOE_%s_%s_%s" % (restraint["dst_segid"], restraint["dst_resid"], restraint["dst_atom"])
            elif per_residue:
                label = "NOE_%s_%s" % (restraint["dst_segid"], restraint["dst_resid"])

            elif restraint["parent"][:4].upper() == "ASSI":
                label = parent_id
            else:
                parent_id = "NOE_%d" % count
                label = parent_id

            try:
                cmd.distance(label, sele1, sele2, quiet=quiet,
                             width=line_width, gap=0, label=0)
            except CmdException:
                print 'FAILED: %s - %s' % (sele1, sele2)
                continue

            cmd.set("dash_color", line_color, label)
        else:
            print line

    cmd.order("NOE*", "yes")

    if not quiet:
        print ' Info: Created distance objects for %d restraints' % (count)

cmd.extend("plot_noe", plot_noe)
