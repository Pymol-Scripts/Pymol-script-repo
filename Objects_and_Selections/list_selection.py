# Using PyMOL commands:
#
list=[]
iterate (name ca),list.append((resi,resn))
print list
 
[('ASP', '1'), ('CYS', '2'), ('ALA', '3'), ('TRP', '4'), ('HIS', '5'), ('LEU',
 '6'), ('GLY', '7'), ('GLU', '8'), ('LEU', '9'), ('VAL', '10'), ('TRP', '11'), 
('CYS', '12'), ('THR', '13')]
 
#!/usr/bin/python
#
# as Python script:
#
from pymol import cmd,stored
stored.list=[]
cmd.iterate("(name ca)","stored.list.append((resi,resn))")
print stored.list
 
[('1', 'ASP'), ('2', 'CYS'), ('3', 'ALA'), ('4', 'TRP'), ('5', 'HIS'), ('6', '
LEU'), ('7', 'GLY'), ('8', 'GLU'), ('9', 'LEU'), ('10', 'VAL'), ('11', 'TRP'), 
('12', 'CYS'), ('13', 'THR')]
 
# Note from Warren: 
#
# The above approach uses a the global pymol variable "stored"
# In recent versions, "cmd.iterate" has been extended to take a dictionary
# as an argument so that you no longer have to use a global variable.
# Avoiding globals helps prevent conflicts between scripts.
#
from pymol import cmd
my_dict = { 'my_list' : [] }
cmd.iterate("(name ca)","my_list.append((resi,resn))",space=my_dict)
print my_dict['my_list']
 
[('1', 'ASP'), ('2', 'CYS'), ('3', 'ALA'), ('4', 'TRP'), ('5', 'HIS'), ('6', '
LEU'), ('7', 'GLY'), ('8', 'GLU'), ('9', 'LEU'), ('10', 'VAL'), ('11', 'TRP'), 
('12', 'CYS'), ('13', 'THR')]

