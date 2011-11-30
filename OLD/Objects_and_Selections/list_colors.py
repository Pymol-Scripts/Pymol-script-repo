#
# This is how to do it from the PyMOL command line or .pml script:
#
iterate all, print color
 
#! /usr/bin/python
#
# and this in a Python script
#
import pymol
pymol.color_list = []
cmd.iterate('all', 'pymol.color_list.append(color)')
print pymol.color_list

