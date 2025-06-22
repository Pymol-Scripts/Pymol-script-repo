'''
See more here: http://www.pymolwiki.org/index.php/Modevectors
'''
from __future__ import print_function
from pymol.cgo import *    # get constants
from math import *
from pymol import cmd


def modevectors(first_obj_frame, last_obj_frame, first_state=1, last_state=1, outname="modevectors", head=1.0, tail=0.3, head_length=1.5, headrgb="1.0,1.0,1.0", tailrgb="1.0,1.0,1.0", cutoff=4.0, skip=0, cut=0.5, atom="CA", stat="show", factor=1.0, notail=0):
    """
    Authors Sean Law & Srinivasa
    Michigan State University
    slaw_(at)_msu_dot_edu

    Editor Sacha Yee

    USAGE

    While in PyMOL

    Parameter                Preset            Type    Description
    first_obj_frame          Undefined         String  Object name of the first structure.  The mode vector will propagate from this structure.  Defined by user.
    last_obj_frame           Undefined         String  Object name of the last structure.  The mode vector (arrow head) will end at this structure.  Defined by user.
    first_state              1                 Integer Defines state of first object
    last_state               1                 Integer Defines state of last object
    outname                  modevectors       String  Name of object to store mode vectors in.
    head                     1.0               Float   Radius for the circular base of the arrow head (cone)
    tail                     0.3               Float   Radius for the cylinder of the arrow tail (cylinder)
    head_length              1.5               Float   Length of the arrow head (from the base of the cone to the tip of cone)
    head_rgb                 1.0,1.0,1.0       String  RGB colour for the arrow head.
    tail_rgb                 1.0,1.0,1.0       String  RGB colour for the arrow tail.
    cutoff                   4.0               Float   Skips mode vectors that do not meet the cutoff distance (in Angstroms).
    skip                     0                 Integer Denotes how many atoms to skip.  No arrows will be created for skipped atoms.
    cut                      0.0               Float   Truncates all arrow tail lengths (without disturbing the arrow head) (in Angstroms).
    atom                     CA                String  Designates the atom to derive mode vectors from.
    stat                     show              String  Keeps track and prints statistics (total modevectors, skipped, cutoff).
    factor                   1.0               Float   Multiplies each mode vector length by a specified factor.
                                                       Values between 0 and 1 will decrease the relative mode vector length.
                                                       Values greater than 1 will increase the relative mode vector length.
    notail                   0                 Integer Hides tails and only uses cones (porcupine plot)
    """

    framefirst = cmd.get_model(first_obj_frame, first_state)
    framelast = cmd.get_model(last_obj_frame, last_state)
    objectname = outname
    factor = float(factor)
    arrow_head_radius = float(head)
    arrow_tail_radius = float(tail)
    arrow_head_length = float(head_length)
    cutoff = float(cutoff)
    skip = int(skip)
    cut = float(cut)
    atomtype = atom.strip('"[]()')
    objectname = objectname.strip('"[]()')

    headrgb = headrgb.strip('" []()')
    tailrgb = tailrgb.strip('" []()')
    hr, hg, hb = list(map(float, headrgb.split(',')))
    tr, tg, tb = list(map(float, tailrgb.split(',')))

    version = cmd.get_version()[1]
    arrow = []
    arrowhead = []
    arrowtail = []
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    z2 = []
    exit_flag = False

##############################################################
#                                                            #
# Define an object called "tail" and store the tail and  a   #
# circular base of the triangle in this object.              #
#                                                            #
##############################################################

    skipcount = 0
    skipcounter = 0
    keepcounter = 0
    atom_lookup = {}
    for atom in framefirst.atom:
        if atom.name == atomtype:
            if skipcount == skip:
                x1.append(atom.coord[0])
                y1.append(atom.coord[1])
                z1.append(atom.coord[2])

                ##########################################
                #                                        #
                # Set atom_lookup for a specific atom    #
                # equal to ONE for the first input set.  #
                # This dictionary will be used as a      #
                # reference for the second set.          #
                #                                        #
                ##########################################

                current_atom = "CHAIN " + atom.chain + " RESID "\
                    + atom.resi + " RESTYPE "\
                    + atom.resn +\
                    " ATMNUM " + str(atom.index)
#				print current_atom
                atom_lookup['current_atom'] = 1

                skipcount = 0
                keepcounter += 1
            else:
                #				print skipcount
                skipcount += 1
                skipcounter += 1

    skipcount = 0
    for atom in framelast.atom:
        if atom.name == atomtype:
            if skipcount == skip:
                x2.append(atom.coord[0])
                y2.append(atom.coord[1])
                z2.append(atom.coord[2])

                #########################################
                #                                       #
                # Get atom information from second set  #
                # and compare with first set.  All      #
                # atoms from this second set MUST be    #
                # found in the first set!  Otherwise,   #
                # the script will exit with an error    #
                # since modevectors can only be created #
                # by calculating values from identical  #
                # sources.                              #
                #                                       #
                #########################################

                current_atom = "CHAIN " + atom.chain + " RESID "\
                    + atom.resi + " RESTYPE "\
                    + atom.resn +\
                    " ATMNUM " + str(atom.index)
#				print current_atom
                if 'current_atom' not in atom_lookup:
                    print("\nError: " + current_atom + " from \""\
                          + last_obj_frame +\
                          " \"is not found in \"" + first_obj_frame + "\".")
                    print("\nPlease check your input and/or selections and try again.")
                    exit_flag = True
                    break

                skipcount = 0
            else:
                skipcount += 1

    if exit_flag == 1:
        ###########################################
        #                                         #
        # Exit script because an atom cannot be   #
        # found in both input files               #
        #                                         #
        ###########################################
        return

    cutoff_counter = 0  # Track number of atoms failing to meet the cutoff

    ###################################################
    #                                                 #
    # Check that the two selections/PDB files contain #
    # the same number of atoms.                       #
    #                                                 #
    ###################################################

    if len(x2) != len(x1):
        print("\nError: \"" + first_obj_frame +\
              "\" and \"" + last_obj_frame +\
              "\" contain different number of residue/atoms.")
        print("\nPlease check your input and/or selections and try again.")
        return
    else:
        # Continue with representing modevectors!
        #########################################
        #                                       #
        # Delete old selection or object if it  #
        # exists so that it can be overwritten  #
        #                                       #
        #########################################
        save_view = cmd.get_view(output=1, quiet=1)
        cmd.delete(objectname)
        cmd.hide(representation="everything", selection=first_obj_frame)
        cmd.hide(representation="everything", selection=last_obj_frame)

    ###################################################
    #                                                 #
    # Begin drawing arrow tails                       #
    #                                                 #
    ###################################################

    arrowtail = []
    for mv in range(len(x1)):
        vectorx = x2[mv] - x1[mv]
        vectory = y2[mv] - y1[mv]
        vectorz = z2[mv] - z1[mv]
        length = sqrt(vectorx ** 2 + vectory ** 2 + vectorz ** 2)
        if length < cutoff:
            cutoff_counter += 1
            continue
        t = 1.0 - (cut / length)
        x2[mv] = x1[mv] + factor * t * vectorx
        y2[mv] = y1[mv] + factor * t * vectory
        z2[mv] = z1[mv] + factor * t * vectorz
        vectorx = x2[mv] - x1[mv]
        vectory = y2[mv] - y1[mv]
        vectorz = z2[mv] - z1[mv]
        length = sqrt(vectorx ** 2 + vectory ** 2 + vectorz ** 2)
        d = arrow_head_length  # Distance from arrow tip to arrow base
        t = 1.0 - (d / length)
        if notail:
            t = 0
        tail = [
            # Tail of cylinder
            CYLINDER, x1[mv], y1[mv], z1[mv]\
            , x1[mv] + (t + 0.01) * vectorx, y1[mv] + (t + 0.01) * vectory, z1[mv] + (t + 0.01) * vectorz\
            , arrow_tail_radius, tr, tg, tb, tr, tg, tb  # Radius and RGB for each cylinder tail
        ]
        if notail == 0:
            arrow.extend(tail)

        x = x1[mv] + t * vectorx
        y = y1[mv] + t * vectory
        z = z1[mv] + t * vectorz
        dx = x2[mv] - x
        dy = y2[mv] - y
        dz = z2[mv] - z
        seg = d / 100
        intfactor = int(factor)
        if version < 1.1:  # Version >= 1.1 has cone primitive
            for i in range(100, 0, -1):  # i=100 is tip of cone
                print(i)
                t1 = seg * i
                t2 = seg * (i + 1)
                radius = arrow_head_radius * (1.0 - i / (100.0))  # Radius of each disc that forms cone
                head = [
                    CYLINDER, x + t2 * dx, y + t2 * dy, z + t2 * dz\
                    , x + t1 * dx, y + t1 * dy, z + t1 * dz\
                    , radius, hr, hg, hb, hr, hg, hb  # Radius and RGB for slice of arrow head
                ]
                arrow.extend(head)
        else:
            head = [
                CONE, x, y, z, x + d * dx, y + d * dy, z + d * dz, arrow_head_radius, 0.0, hr, hg, hb, hr, hg, hb, 1.0, 1.0]
            arrow.extend(head)

##############################################################
#                                                            #
# Load the entire object into PyMOL                          #
#                                                            #
# Print statistics if requested by user                      #
#                                                            #
##############################################################

    if stat == "show":
        natoms = skipcounter + keepcounter
        print("\nTotal number of atoms = " + str(natoms))
        print("Atoms skipped = " + str(skipcounter))
        if keepcounter - cutoff_counter > 0:
            print("Atoms counted = " + str(keepcounter - cutoff_counter) + " (see PyMOL object \"" + objectname + "\")")
        else:
            print("Atoms counted = " + str(keepcounter - cutoff_counter) + " (Empty CGO object not loaded)")
        print("Atoms cutoff  = " + str(cutoff_counter))  # Note that cutoff occurs AFTER skipping!
    if keepcounter - cutoff_counter > 0:
        cmd.delete(objectname)
        cmd.load_cgo(arrow, objectname)  # Ray tracing an empty object will cause a segmentation fault.  No arrows = Do not display in PyMOL!!!
    cmd.show(representation="cartoon", selection=first_obj_frame)
    if (first_obj_frame != last_obj_frame):
        cmd.show(representation="cartoon", selection=last_obj_frame)
        cmd.hide(representation="cartoon", selection=last_obj_frame)
    cmd.bg_color(color="white")
    cmd.set_view(save_view)
    return

cmd.extend("modevectors", modevectors)
