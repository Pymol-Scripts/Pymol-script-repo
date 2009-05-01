# This script writes the distance from 
# atom mol1///25/ha to atom mol1///26/ha
# out to the file "dist.txt"
# Simply change your selections to see different distances.
 
# import PyMOL's command namespace
from pymol import cmd
 
# open dist.txt for writing
f=open('dist.txt','w')
 
# calculate the distance and store it in dst
dst=cmd.distance('tmp','mol1///25/ha','mol1///26/ha')
 
# write the formatted value of the distance (dst)
# to the output file
f.write("%8.3f\n"%dst)
 
# close the output file.
f.close()
