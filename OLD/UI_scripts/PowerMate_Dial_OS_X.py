from pymol import cmd
 
# Define aliases for mapping in [x,y,z] rotations and translations into a single Powermate
# dial.  Toggling between the three is possible if you then assign these to special keys.
 
# Functions for x,y,z rotations and translations using Powermate dial
# Program F1 and F2 for Rotate Right and Rotate Left
# Program F3 and F4 for Click & Rotate Right and Click & Rotate Left
# Program F5 for  Click  (to toggle between dialsets)
 
# dialset = 2
 
def dialx(): \
    global dialset \
    dialset = 1 \
    cmd.set_key ('F1', cmd.turn,('x',-2.0)) \
    cmd.set_key ('F2', cmd.turn,('x',2.0)) \
    cmd.set_key ('F3', cmd.move,('x',-0.5)) \
    cmd.set_key ('F4', cmd.move,('x',0.5)) \
    print "dialset ", dialset, " [ X ]\n" \
    return dialset
 
def dialy(): \
    global dialset \
    dialset = 2 \
    cmd.set_key ('F1', cmd.turn,('y',-2.0)) \
    cmd.set_key ('F2', cmd.turn,('y',2.0)) \
    cmd.set_key ('F3', cmd.move,('y',-0.5)) \
    cmd.set_key ('F4', cmd.move,('y',0.5)) \
    print "dialset ", dialset, " [ Y ]\n" \
    return dialset
 
 
def dialz(): \
    global dialset \
    dialset = 3 \
    cmd.set_key ('F1', cmd.turn,('z',-2.0)) \
    cmd.set_key ('F2', cmd.turn,('z',2.0)) \
    cmd.set_key ('F3', cmd.move,('z',-0.5)) \
    cmd.set_key ('F4', cmd.move,('z',0.5)) \
    print "dialset ", dialset, " [ Z ]\n" \
    return dialset
 
def toggle_dial(): \
    if dialset == 1 : \
        print "Changing to y" \
        dialy() \
    elif dialset == 2 : \
        print "Changing to z" \
        dialz() \
    elif dialset == 3 : \
        print "Changing to x" \
        dialx() \
    else: print "Dial assignment isn't working"
 
 
cmd.set_key ('F5', toggle_dial)
 
# Start default dial state for rotate y  (arbitrary choice)
 
dialy()
