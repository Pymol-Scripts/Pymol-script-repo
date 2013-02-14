'''
http://pymolwiki.org/index.php/plot_noe

(c) 2011-09-23 by Justin L Lorieau
(c) 2012-10-31 by Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def plot_noe(filename, selection='', line_color='gray20', line_width='1.0', single=0, quiet=1):
    """
DESCRIPTION

    A function for plotting XPLOR NOE restraints on a structure

ARGUMENTS

    filename = string: The filename of the NOE retraint file in XPLOR NIH format.

    selection = string: atom selection {default: all}

    line_color = string: The color for the NOE lines. {default: black}

    line_width = float: The thickness of the NOE lines. {default: 1.0}

NOE Restraint Format Example

    assign (residue 5 and name HB#) (residue 21 and name HA) 3.0 0.7 0.7

EXAMPLE

    PyMOL> plot_noe noe_short.tbl
    """
    from pymol.parsing import split

    single, quiet = int(single), int(quiet)

    count = 0

    for line in open(filename):
        line = line.replace('(', ' (').replace(')', ') ').replace('#', '*')

        a = filter(None, split(line, ' \t'))
        if len(a) < 4 or not a[0].startswith('assi'):
            continue

        try:
            float(a[3])
        except ValueError:
            continue

        sele1, sele2 = a[1:3]

        if selection:
            sele1 = '(%s) and (%s)' % (selection, sele1)
            sele2 = '(%s) and (%s)' % (selection, sele2)

        label = single and "NOE" or ("NOE_%d" % (count + 1))

        try:
            cmd.distance(label, sele1, sele2, quiet=quiet,
                    width=line_width, gap=0, label=0)
        except CmdException:
            print 'FAILED: %s - %s' % (sele1, sele2)
            continue

        cmd.set("dash_color", line_color, label)
        count += 1

    if not quiet:
        print ' Info: Created distance objects for %d restraints' % (count)

cmd.extend("plot_noe", plot_noe)

