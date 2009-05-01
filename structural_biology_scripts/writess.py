import os;
import os.path;
import glob;
import string;
 
from pymol import cmd
from pymol import stored
from pymol import selector
 
files = glob.glob("/tmp/thy_model/*.pdb")
 
for file in files:
        pdbName = string.split(os.path.basename(file), ".")[0]
        cmd.load(file, pdbName)
        outFile = open(pdbName + '.ss', 'wb')
        stored.ss = ""
        cmd.iterate( '(n. CA)', 'stored.ss = stored.ss + ("%1s"%ss)')
        for c in stored.ss:
                if c  == " ":
                        outFile.write('.')
                else:
                        outFile.write(c)
        cmd.delete(pdbName)
        outFile.close()

