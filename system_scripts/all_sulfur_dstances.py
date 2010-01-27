from pymol import cmd
from glob import glob
 
for file in glob("*.pdb"):
    print file
    cmd.load(file,'prot')
    for a in cmd.index("CYS/SG"):
        for b in cmd.index("CYS/SG"):
            if a[1]<b[1]:
                cmd.select("s1","%s`%d"%a)
                cmd.select("s2","%s`%d"%b)
                if cmd.select("(s1|s2) and not ?skip"):
                    cmd.iterate("s1|s2","print ' ',chain,resn,resi,name")
                    print '   ',round(cmd.dist("tmp","s1","s2"),3)
                    cmd.select("skip","s1|s2|?skip")
    cmd.delete("all")
