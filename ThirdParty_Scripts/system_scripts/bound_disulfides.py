from pymol import cmd
from glob import glob
 
for file in glob("*.pdb"):
    print file
    cmd.load(file,'prot')
    for a in cmd.index("elem s and bound_to elem s"):
        if cmd.select("s1","%s`%d"%a) and \
           cmd.select("s2","elem s and bound_to %s`%d"%a):
            if cmd.select("(s1|s2) and not ?skip"):
                cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
                print '   ',round(cmd.dist("tmp","s1","s2"),3)
                cmd.select("skip","s1|s2|?skip")
    cmd.delete("all")

