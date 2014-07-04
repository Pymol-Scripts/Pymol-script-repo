'''
http://www.pymolwiki.org/index.php/viol_noes
 
(c) August 2010 by Mateusz Maciejewski
matt (at) mattmaciejewski . com

License: MIT

'''


from pymol.cgo import *
from pymol import cmd


def viol_noes(molecule, viol_file, viol_class=None, quiet=1):
    """
DESCRIPTION

    Visualize Xplor-NIH NOE violations.

ARGUMENTS

    molecule = string: molecule on which to show the violations.

    viol_file = string: Xplor-NIH .viol file that contains the violations to be visualized.

    viol_class = string: NOE class in .viol file to show {default: None (means all NOE classes)}.


EXAMPLE

    PyMOL> run viol_noes.py
    PyMOL> viol_noes molecule, ./molecule.pdb.viol

NOTES

    The NOE violations will be shown as distances between the relevant residues/atoms and colored according to the severity of violation (the closer to the blue end of the spectrum, the more severe the violation; to closer to the red end, the less severe the violation).
    """

    import sys
    import re

    status = None
    classMk = None

    workfile = open(viol_file, 'r')

    for line in workfile.readlines():

        line = line.strip()

        if line.find('Violated NOE') != -1:
            classMk = line.split()[6]
            status = "new_class"

        elif not line.strip():
            pass

        elif line == "number of restraints":
            status = "end_class"

        elif status == "new_class" and line.find("----") != -1:
            status = "in_class"

        elif status == "in_class" and (viol_class == None or classMk == viol_class):

            try:

                new_id = float(line.split()[0])
                alt = 1  # in case there is a couple of alternatives this labels them
                line = line.replace('(', '').replace(')', '').split()
                member1 = (line[1], line[3])
                member2 = (line[4], line[6])
                delta = line[9]

                cmd.distance("viol_%d_%s_%s_%d" % (new_id, member1[0], member2[0], alt), " %s & i. %s & n. %s" % (molecule, member1[0], member1[1]), "%s & i. %s & n. %s" % (molecule, member2[0], member2[1]), quiet=int(quiet))
                cmd.color("o%s" % (int((100) * float(delta))), "viol_%d_%s_%s_%d" % (new_id, member1[0], member2[0], alt))

            except ValueError:

                alt += 1

                if line.split()[1] != ")":
                    if line.split()[-2] == "(" and line.split()[-1] == ")":
                        line = line.replace('(', '').replace(')', '').split()
                        member1 = (line[0], line[2])

                elif line.split()[1] == ")":
                    line = line.replace('(', '').replace(')', '').split()
                    member2 = (line[-3], line[-1])

                else:
                    line = line.replace('(', '').replace(')', '').split()
                    member1 = (line[0], line[2])
                    member2 = (line[-3], line[-1])

                cmd.distance("viol_%d_%s_%s_%d" % (new_id, member1[0], member2[0], alt), " %s & i. %s & n. %s" % (molecule, member1[0], member1[1]), "%s & i. %s & n. %s" % (molecule, member2[0], member2[1]), quiet=int(quiet))
                cmd.color("o%s" % (int((100) * float(delta))), "viol_%d_%s_%s_%d" % (new_id, member1[0], member2[0], alt))

    cmd.hide("labels", "all")
    cmd.set("dash_radius", 0.03)


cmd.extend("viol_noes", viol_noes)
