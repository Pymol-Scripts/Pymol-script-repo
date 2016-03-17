'''
This PyMOL module is described at: http://www.pymolwiki.org/Cgo_grid
################################################################################

Author : Andreas Warnecke
email: 4ndreas.warneck3@gmail.com
Date: June 2014
License: pending...
Citation: pending...
Version: 1.0

Module contributed by Andreas Warnecke
(andreas.warnecke@ki.se, 4ndreas.warneck3@gmail.com)

Feel free to contact me in case of feedback (suggestions/comments) or questions.

cgo_grid renders a mesh-like grid with overlapping sine waves (cf. examples).

################################################################################
'''
from __future__ import print_function
from pymol import cmd
from pymol import stored
from pymol.cgo import *
from chempy import cpv
import random
import math

def cgo_grid(
pos1=[0,0,0],
pos2=[1,0,0],
pos3=[0,0,1],
length_x=30,
length_z='',
npoints_x='',
npoints_z='',
nwaves_x=2,
nwaves_z='',
offset_x=0,
offset_z='',
gain_x=1,
gain_z='',
thickness=2.0,
color='',
nstates=60,
startframe=1,
endframe=1,
mode=0,
view=0,
name='',
quiet=1):
    '''
DESCRIPTION

    Generates an animated flowing mesh object using the points provided
    or the current view. The shape is affected substantially by the arguments!

USEAGE

    cgo_grid [ pos1 [, pos2 [, pos3 [, length_x [, length_z
             [, npoints_x [, npoints_z [, nwaves_x [, nwaves_z
             [, offset_x [, offset_z [, gain_x [, gain_z [, thickness
             [, color [, nstates [, startframe [, endframe [, mode
             [, view [, name [, quiet ]]]]]]]]]]]]]]]]]]]]]]

EXAMPLE

    cgo_grid view=1

ARGUMENTS

    pos1 = single atom selection (='pk1') or list of 3 floats {default: [0,0,0]}

    pos2 = single atom selection (='pk2') or list of 3 floats {default: [1,0,0]}

    pos3 = single atom selection (='pk3') or list of 3 floats {default: [0,0,1]}

    --> the plane is defined by pos1 (origin) and vectors to pos2 and pos3, respectively

    length_x = <float>: length of membrane {default: 30}
    length_z = <float>: length of membrane {default: ''} # same as length_x

    npoints_x = <int>: number of points(lines) along x-direction
                {default: ''} #will be set to give a ~1 unit grid
    npoints_z = <int>: number of points(lines) along z-direction
                {default: ''} #will be set to give a ~1 unit grid
                {minimum: 1 # automatic}

    nwaves_x =   <float>: number of complete sin waves along object x-axis
                 {default: 2}
    nwaves_z =  <float>: number of complete sin waves along object z-axis
                {default: ''} # same as nwaves_x
                define separately to adjust number of waves in each direction



    offset_x = <float> phase delay (in degrees) of sin wave in x-axis
             can be set to affect shape and starting amplitude {default: 0}
    offset_z = <float> phase delay (in degrees) of sin wave in z-axis
             can be set to affect shape and starting amplitude
             {default: ''} # same as  offset_x
    offset_x and offset_z can be used together to phase
    otherwise identical objects

    gain_x = <float>: multiplication factor for y-amplitude for x-direction
             {default: 1}
    gain_z = <float>: multiplication factor for y-amplitude for z-direction
             {default: ''} #=gain_x

    thickness = <float>: line thickness {default: 2}

    color = color name <string> (e.g. 'skyblue') OR
            rgb-value list of 3 floats (e.g. [1.0,1.0,1.0]) OR
            {default: ''} // opposite of background
            input illegal values for random coloring

    nstates =  <int>: number of states; {default: 60}
               this setting will define how many states
               the object will have (per wave) and how fluent and fast the
               animation will be.
               Higher values will promote 'fluent' transitions,
               but decrease flow speed.
                   Note: Frame animation cycles thought the states one at a time
                   and needs to be set accordingly. Can also be used to phase
                   otherwise identical objects.
               Set to 1 for static object {automatic minimum}

    startframe: specify starting frame <int> or set (='') to use current frame
                set to 'append' to extend movie from the last frame {default: 1}
      endframe: specify end frame <int> or set (='') to use last frame
                if 'append' is used for startframe,
                endframe becomes the number of frames to be appended instead
                {default: 1}
                Note: if start- and endframe are the same, movie animation will
                be skipped, the object will be loaded and can be used afterwards

    mode: defines positioning {default: 0}:
    0: pos1 is center
    1: pos1 is corner

    view {default: 0}:
    '0': off/ uses provided points to create CGO
    '1': overrides atom selections and uses current orienatation for positioning
         - pos1 = origin/center
         - pos2 = origin +1 in camera y
         - pos3 = origin +1 in camera z

    name: <string> name of cgo object {default: ''} / automatic

    quiet: <boolean> toggles output

    '''
    ########## BEGIN OF FUNCTION CODE ##########
    def get_coord(v):
        if not isinstance(v, str):
            try:
                return v[:3]
            except:
                return False
        if v.startswith('['):
            return cmd.safe_list_eval(v)[:3]
        try:
            if cmd.count_atoms(v)==1:
                # atom coordinates
                return cmd.get_atom_coords(v)
            else:
                # more than one atom --> use "center"
                # alt check!
                if cmd.count_atoms('(alt *) and not (alt "")')!=0:
                    print("cgo_grid: warning! alternative coordinates found for origin, using center!")
                view_temp=cmd.get_view()
                cmd.zoom(v)
                v=cmd.get_position()
                cmd.set_view(view_temp)
                return v
        except:
            return False

    def eval_color(v):
        try:
            if not v:
                v=eval(cmd.get('bg_rgb'))
                v=list(map(sum, list(zip(v,[-1,-1,-1]))))
                v=list(map(abs, v))
                if v[0]==v[1]==v[2]==0.5: # grey
                    v=[0,0,0]
                return v
            if isinstance(v, list):
                return v[0:3]
            if not isinstance(v, str):
                return v[0:3]
            if v.startswith('['):
                return cmd.safe_list_eval(v)[0:3]
            return list(cmd.get_color_tuple(v))
        except:
            return [random.random(),random.random(),random.random()]
    cmd.extend("eval_color", eval_color)

    color=eval_color(color)

    try:
        mode=int(mode)
    except:
        raise Exception("Input error in Mode")
    if mode<0 or mode>1:
        raise Exception("Mode out of range!")

    try:
        nstates=int(nstates)
        if nstates<1:
            nstates=1
            print("NB! nstates set to 1 (automatic minimum)")
        length_x=float(length_x)
        if length_z=='':
            length_z=length_x
        else:
            length_z=float(length_z)
        if npoints_x=='':
            npoints_x=int(length_x)+1
        else:
            npoints_x=int(npoints_x)
        if npoints_x<1:
            npoints_x=1
            print("NB! npoints_x set to 1 (automatic minimum)")
        if npoints_z =='':
            npoints_z=int(length_z)+1
        else:
            npoints_z=int(npoints_z)
        if npoints_z<1:
            npoints_z=1
            print("NB! npoints_x set to 1 (automatic minimum)")

        nwaves_x=abs(float(nwaves_x))
        if nwaves_z=='':
            nwaves_z=nwaves_x
        else:
            nwaves_z=abs(float(nwaves_z))
        offset_x=float(offset_x)*math.pi/180
        if offset_z=='':
            offset_z=offset_x
        else:
            offset_z=float(offset_z)*math.pi/180
        thickness=float(thickness)
        gain_x=float(gain_x)
        if gain_z=='':
            gain_z=gain_x
        else:
            gain_z=float(gain_z)
        if not name:
            name = cmd.get_unused_name('membrane')
        else:
            name = str(name)

        if int(quiet):
            quiet=True
        else:
            quiet=False
        if int(view):
            view=True
        else:
            view=False
    except:
        raise Exception("Input error in parameters!")


    #prevent auto zooming on object
    temp_auto_zoom=cmd.get('auto_zoom')
    cmd.set('auto_zoom', '0')

    if int(view):
        xyz1=cmd.get_position()
        tempname = cmd.get_unused_name('temp')
        ori_ax=[[0,0,0],[10,0,0],[0,0,10]]
        for a in range (0,len(ori_ax)):
            cmd.pseudoatom(tempname, resi=''+str(a+1)+'', pos=xyz1)
            cmd.translate(ori_ax[a],
            selection=''+tempname+' and resi '+str(a+1)+'', camera='1')
            ori_ax[a]=cmd.get_atom_coords(''+tempname+' and resi '+str(a+1)+'')
        cmd.delete(tempname)
        xyz1=ori_ax[0]
        xyz2=ori_ax[1]
        xyz3=ori_ax[2]
    else:
        xyz1 = get_coord(pos1)
        xyz2 = get_coord(pos2)
        xyz3 = get_coord(pos3)

    if (not startframe):
        startframe=cmd.get('frame')

    if (not endframe):
        endframe=cmd.count_frames()
    if endframe==0: endframe=1

    if (startframe=='append'):
        startframe=cmd.count_frames()+1
        try:
            endframe=int(endframe)
            cmd.madd('1 x'+str(endframe))
            endframe=cmd.count_frames()
        except ValueError:
            raise Exception("Input error: Value for 'endframe' is not integer!")

    try:
        startframe=int(startframe)
        endframe=int(endframe)
        endframe/startframe
        startframe/endframe
    except ValueError:
        raise Exception("Input error: Failed to convert to integer!")
    except ZeroDivisionError:
        raise Exception("Error: unexpected zero value!")
    except:
        raise Exception("Unexpected error!")

    if (nstates==1):
        if not quiet: print("Creating one state object!")

    if startframe > endframe:
        startframe, endframe = endframe, startframe
        if not quiet: print("Inverted start and end frames!")


    ########## BEGIN OF FUNCTIONAL SCRIPT ##########

    #normalize and get orthogonal vector

    # define vectors from points
    xyz2 = cpv.sub(xyz2, xyz1)
    xyz3 = cpv.sub(xyz3, xyz1)

    #NB! cpv.get_system2 outputs normalized vectors [x,y,z]
    xyz4 = cpv.get_system2(xyz2,xyz3)
    xyz2 = xyz4[0]
    xyz3 = xyz4[1]
    for x in range(0,3):
        for z in range(0,3):
            if x==z:
                continue
            if xyz4[x]==xyz4[z]:
                raise Exception("Illegal vector settings!")
    xyz4 = cpv.negate(xyz4[2]) #overwrites original

    # transform origin to corner
    if mode==0:
        if npoints_x>1:
            xyz1 = cpv.sub(xyz1, cpv.scale(xyz2,length_x/2))
        if npoints_z>1:
            xyz1 = cpv.sub(xyz1, cpv.scale(xyz3,length_z/2))

    #defines array lines
    nlines=max([npoints_x, npoints_z])
    # in case only one line max

    # create an empty array for xyz entries
    # this may contain more values than are actually drawn later,
    # but they are needed to draw lines in each direction
    grid_xyz = []
    for x in range(0,nlines):
        grid_xyz.append([0.0,0.0,0.0]*nlines)

    # grid distance and steps
    # prevent zero divisions (lines=1) and enable calculations if lines=0
    if (not (npoints_x-1<2)):
        gap_length_x = length_x/(npoints_x-1)
        step_line_x = 2*math.pi/(npoints_x-1)
    else:
        gap_length_x=length_x
        step_line_x=2*math.pi

    if (not (npoints_z-1<2)):
        gap_length_z = length_z/(npoints_z-1)
        step_line_z = 2*math.pi/(npoints_z-1)
    else:
        gap_length_z=length_z
        step_line_z=2*math.pi

    # calculate steps
    if nstates==1:
        step_state=0
    else:
        step_state = 2*math.pi/(nstates-1)

    ########## BEGIN STATE ITERATION ##########
    # create a n-state object in PyMol

    for a in range(0,nstates):
        # Reset object
        obj = []
        #assign color
        obj.extend( [ COLOR, color[0], color[1], color[2] ] )
        #set width
        obj.extend( [ LINEWIDTH, thickness ] )

        # Calculate xyz-coordinates for each line

        for x in range(0,nlines):

            for z in range(0,nlines):

                # update grid position in x-direction
                xyztemp=cpv.add(xyz1,cpv.scale(xyz2,gap_length_x*x))

                # update grid position in z-direction
                xyztemp=cpv.add(xyztemp,cpv.scale(xyz3,gap_length_z*z))

                # calculate amplitude for y-direction and update grid position
                y_amp=(\
                      gain_x*math.sin(offset_x+nwaves_x*((a*step_state)+(x*step_line_x)))/2+\
                      gain_z*math.sin(offset_z+nwaves_z*((a*step_state)+(z*step_line_z)))/2\
                      )
                xyztemp=cpv.add(xyztemp,cpv.scale(xyz4,y_amp))
                grid_xyz[x][z]=xyztemp

        #Now the coordinates for this state are defined!

        #Now the coordinates are read separately:

        # allow to run the loops as often as required
        #if npoints_x==0:npoints_x=npoints_z

        #lines along z in x direction
        for z in range(0,npoints_z):
            obj.extend( [ BEGIN, LINE_STRIP ] )
            for x in range(0,npoints_x):
                obj.extend( [ VERTEX, grid_xyz[x][z][0], grid_xyz[x][z][1], grid_xyz[x][z][2] ] )
            obj.append( END )

        #lines along x in z direction
        for x in range(0,npoints_x):
            obj.extend( [ BEGIN, LINE_STRIP ] )
            for z in range(0,npoints_z):
                obj.extend( [ VERTEX, grid_xyz[x][z][0], grid_xyz[x][z][1], grid_xyz[x][z][2] ] )
            obj.append( END )

        # Load state into PyMOL object:
        cmd.load_cgo(obj,name,a+1)
    # All states of object loaded!

    #reset auto zooming to previous value
    cmd.set('auto_zoom', temp_auto_zoom)

    # animate object using frames instead of states
    if (not endframe==startframe):
        framecount=0
        countvar=1

        for frame in range(startframe, endframe + 1):
            #increase count
            framecount=framecount+countvar

            # set state in frame
            cmd.mappend(frame,
            "/cmd.set('state', %s, %s)" % (repr(framecount), repr(name)))

            # Looping
            if framecount==nstates:
                if ((int(nwaves_x)!=nwaves_x) or (int(nwaves_z)!=nwaves_z)):
                    #if not complete sinus wave
                    #--> reverse wave in ongoing animation
                    countvar=-1
                else:
                    #wave is complete --> repeat
                    framecount=0
            # count up from first state
            if framecount==1: countvar=1
        if not quiet: print("object loaded and animated with frames!")
    else:
        if not quiet: print("object loaded!")

    #OUTPUT
    if not quiet:
        print("Grid variables for:",name)
        print("corner:", xyz1)
        print("vector 1:", xyz2)
        print("vector 2:", xyz3)

        print("length_x:",length_x)
        print("length_z:",length_z)
        print("npoints_x:", npoints_x)
        print("npoints_z:", npoints_z)

        print("nwaves_x:", nwaves_x)
        print("nwaves_z:", nwaves_z)

        print("offset_x:",offset_x)
        print("offset_z:",offset_z)

        print("gain_x:",gain_x)
        print("gain_z:",gain_z)

        print("thickness:",thickness)

        print("states", nstates)
        if (not endframe==startframe):
            print("frames: start:",startframe,"end:",endframe)

    return grid_xyz

cmd.extend("cgo_grid", cgo_grid)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
