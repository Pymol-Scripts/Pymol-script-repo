"""
http://pymolwiki.org/index.php/plot_noe

(c) 2011-09-23 by Justin L Lorieau
(c) 2012-10-31 by Thomas Holder

License: BSD-2-Clause
"""

from __future__ import print_function

from pymol import cmd, CmdException
import re
import shlex


def plot_noe(filename, line_color=None, line_width='1.0', advanced_coloring=0, single=0, quiet=1, aria=0, per_atom=0,
             per_residue=1):
    """
DESCRIPTION

    A function for plotting XPLOR NOE restraints on a structure

ARGUMENTS

    filename = string: The filename of the NOE retraint file in XPLOR NIH format.

    line_color = string: The color for the NOE lines. {default: yellow}

    line_width = float: The thickness of the NOE lines. {default: 1.0}

    advanced_coloring = color restraints by distance.

    single = string: create a single object for all restraints.

    aria = integer: Name NOEs after Aria IDs.

    per_atom: Group NOEs on atom basis

    per_residue: Group NOEs on residue basis (default)

NOE Restraint Format Example

    assign (residue 5 and name HB#) (residue 21 and name HA) 3.0 0.7 0.7

EXAMPLE

    PyMOL> plot_noe noe_short.tbl
    """

    single, quiet, aria, per_residue, per_atom, advanced_coloring = \
        int(single), int(quiet), int(aria), int(per_residue), int(per_atom), int(advanced_coloring)

    count = 0

    if advanced_coloring == 1:
        cmd.set("group_auto_mode", 2)

    rgx_assi = re.compile("\s*[asignASIGN]+\s+.*")
    rgx_or = re.compile("\s*[orOR]+\s+.*")
    rgx_bracket = re.compile("(\(|\))")
    rgx_tick = re.compile("(\w)\'")

    restraint_line = None
    restraint_block = []
    restraints_blocks = []

    for line in open(filename):
        if rgx_assi.match(line):
            if restraint_line:
                restraint_block.append(rgx_tick.sub("\g<1>^", rgx_bracket.sub(" ", restraint_line)))
                restraints_blocks.append(restraint_block)
            restraint_block = []
            restraint_line = line
        elif rgx_or.match(line):
            restraint_block.append(rgx_tick.sub("\g<1>^", rgx_bracket.sub(" ", restraint_line)))
            restraint_line = line
        else:
            restraint_line += line

    restraint_block.append(rgx_tick.sub("\g<1>^", rgx_bracket.sub(" ", restraint_line)))
    restraints_blocks.append(restraint_block)

    for restraint_block in restraints_blocks:
        restraints = {"connections": []}

        try:
            shlex_split = shlex.split(restraint_block[0])
        except ValueError:
            print('Cannot process "%s"' % restraint_block[0])
        try:
            restraints["distance"] = float(shlex_split[11].strip())
        except IndexError:
            print('Cannot process "%s"' % "".join(restraint_block))
            continue
        except ValueError:
            try:
                restraints["distance"] = float(shlex_split[17].strip())
            except ValueError:
                print("Failed to extract distance, setting to 1A")
                restraints["distance"] = 1

        if len(shlex_split) > 20:
            restraints["ID"] = shlex_split[26][3:-1].strip()
        else:
            restraints["ID"] = str(count + 1)

        for restraint_line in restraint_block:
            distance = {}
            try:
                shlex_split = shlex.split(restraint_line)

                if len(shlex_split) > 6:
                    if shlex_split[1].strip().lower() == shlex_split[9].strip().lower() == "segid":
                        distance["src_segid"] = shlex_split[2].strip()
                        distance["src_resid"] = shlex_split[5].strip()
                        distance["src_atom"] = shlex_split[8].strip().replace("^", "'")
                        distance["dst_segid"] = shlex_split[10].strip()
                        distance["dst_resid"] = shlex_split[13].strip()
                        distance["dst_atom"] = shlex_split[16].strip().replace("^", "'")
                    else:
                        distance["src_segid"] = "A"
                        distance["src_resid"] = shlex_split[2].strip()
                        distance["src_atom"] = shlex_split[5].strip().replace("^", "'")
                        distance["dst_segid"] = "A"
                        distance["dst_resid"] = shlex_split[10].strip()
                        distance["dst_atom"] = shlex_split[13].strip().replace("^", "'")

                else:
                    print('Cannot process "%s"' % restraint_line)
                    continue
            except IndexError:
                print('Cannot process "%s"' % restraint_line)
                continue
            except ValueError:
                print('Cannot process "%s"' % restraint_line)
                continue

            restraints["connections"].append(distance)

        count += _draw_restraint(
            restraints,
            line_color,
            line_width,
            single,
            quiet,
            per_residue,
            per_atom,
            advanced_coloring
        )

    if not quiet:
        print(' Info: Created distance objects for %d restraints' % count)

    cmd.order("NOE*", "yes")


def _draw_restraint(restraints, line_color, line_width, single, quiet, per_residue, per_atom, advanced_coloring):
    cnt = 0
    for restraint in restraints["connections"]:
        cnt += 1
        sele1 = 'segid %s and resi %s and name %s' % (
            restraint["src_segid"], restraint["src_resid"], restraint["src_atom"].replace("#", "*")
        )
        sele2 = 'segid %s and resi %s and name %s' % (
            restraint["dst_segid"], restraint["dst_resid"], restraint["dst_atom"].replace("#", "*")
        )

        if line_color is None:
            if advanced_coloring == 1:
                if restraints["distance"] < 2 or restraints["distance"] > 6:
                    line_color = "red"
                elif 5 < restraints["distance"] < 6:
                    line_color = "orange"
                else:
                    line_color = "splitpea"
            else:
                line_color = "yellow"

        if per_atom or per_residue:
            for x, y, z in [
                [restraint["dst_segid"], restraint["dst_resid"], restraint["dst_atom"]],
                [restraint["src_segid"], restraint["src_resid"], restraint["src_atom"]]
            ]:
                if per_atom:
                    label = "NOE_%s_%s_%s" % (x, y, z)
                elif per_residue:
                    label = "NOE_%s_%s" % (x, y)
                if advanced_coloring == 1:
                    label += ".%s_%i" % (restraints["ID"], cnt)

                try:
                    cmd.distance(label, sele1, sele2, quiet=quiet,
                                 width=line_width, gap=0, label=0)
                    _color_restraint(label, line_color)
                except CmdException:
                    print('FAILED: %s - %s' % (sele1, sele2))
                    print(restraints["ID"])
                    return 0

        else:
            if single:
                label = "NOE"

            else:
                label = "NOE_" + restraints["ID"]

            try:
                cmd.distance(label, sele1, sele2, quiet=quiet,
                             width=line_width, gap=0, label=0)

                _color_restraint(label, line_color)
            except CmdException:
                print('FAILED: %s - %s' % (sele1, sele2))
                print(restraints["ID"])
                return 0

    return 1


def _color_restraint(object_name, color="yellow"):
    cmd.set("dash_color", color, object_name)


cmd.extend("plot_noe", plot_noe)
