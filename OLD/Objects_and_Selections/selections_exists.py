#
# The function below will return true if the selection is defined.
#
from pymol import cmd
from types import *
 
def sele_exists(sele):
    sess = cmd.get_session()
    for i in sess["names"]:
        if type(i) is ListType:
            if sele==i[0]:
                return 1
    return 0
 
# Note from Warren:
#
# simpler, faster:
 
def sele_exists(sele):
   return sele in cmd.get_names("selections")

