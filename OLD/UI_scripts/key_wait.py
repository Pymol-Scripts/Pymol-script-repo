#use "spawn spawn_demo.py, local" to invoke this python script from
within pymol
wait=""
i=0
while i < 12 and wait!="x":
  cmd.turn("y", 30)
  print "Press enter key to continue or x + enter to terminate"
  wait=raw_input()
  i=i+1
print "Done"
