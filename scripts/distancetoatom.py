'''
distancetoatom.py
Described at: http://www.pymolwiki.org/Distancetoatom
##################################################

Prints all distanced between the specified atom/coordinate/center
and all atoms within cutoff distance that are part of the selection.
All coordinates and distances can be saved in a csv-style text file report
and can be appended to a (custom) atom property, if defined.


##################################################
Author : Andreas Warnecke, Jared Sampson
email: 4ndreas.warneck3@gmail.com, Jared.Sampson@nyumc.org
Date: June 2014
License: BSD2
Version: 1.0
##################################################
VERSION NOTES:
    1.0    2014    First release
'''

from __future__ import print_function
import sys
from pymol import cmd
from pymol import stored
from chempy import cpv
from pymol import cmd

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
                print("distancetoatom: warning! alternative coordinates found for origin, using center!")
            view_temp=cmd.get_view()
            cmd.zoom(v)
            v=cmd.get_position()
            cmd.set_view(view_temp)
            return v
    except:
        return False


def distancetoatom(
origin='pk1',
cutoff=10,
filename=None,
selection='all',
state=0,
property_name='p.dist',
coordinates=0,
decimals=3,
sort=1,
quiet=1
):
    '''
DESCRIPTION

    distancetoatom.py
    Described at: http://www.pymolwiki.org/Distancetoatom

    Prints all distanced between the specified atom/coordinate/center
    and all atoms within cutoff distance that are part of the selection.
    All coordinates and distances can be saved in a csv-style text file report
    and can be appended to a (custom) atom property, if defined.

USAGE

    distancetoatom [ origin [, cutoff [, filename [, selection
    [, state [, property_name [, coordinates [, decimals [, sort
    [, quiet ]]]]]]]]]]

ARGUMENTS

    NAME        TYPE    FUNCTION
    origin:     <list>  defines the coordinates for the origin and can be:
                <str>   1. a list with coordinates [x,y,z]
                        2. a single atom selection string {default='pk1'}
                        3. a multi-atom selection string (center will be used)
    cutoff      <float> sets the maximum distance {default: 10}
    filename    <str>   filename for optional output report. {default=None}
                        set to e.g. 'report.txt' to create a report
                        (omit or set to '', None, 0 or False to disable)
    selection   <str>   can be used to define/limit the measurment to specific
                        sub-selections {default='all'}
    state       <int>   object state, {default=0} # = current
    property_name <str> the distance will be stored in this property {p.dist}
                        set "" to disable
    coordinates <int>   toggle whether atom coordinated will be reported {0}
    decimals    <int>   decimals for coordinates and distance:
                        default = 3 # = max. PDB resolution
    sort        <int>   Sorting by distance?
                         1: ascending (default)
                         0: no sorting (by names)
                        -1: descending
    quiet       <bool>  toggle verbosity
    '''
    # keyword check
    try:
        selection = '(%s)'%selection
        ori=get_coord(origin)
        if not ori:
            print("distancetoatom: aborting - check input for 'origin'!")
            return False
        cutoff = abs(float(cutoff))
        filename = str(filename)
        state = abs(int(state))
        if (not int(state)):
            state=cmd.get_state()
        cmd.set('state', state) # distance by state
        property_name = str(property_name)
        decimals = abs(int(decimals))
        sort = int(sort)
        coordinates=bool(int(coordinates))
        quiet=bool(int(quiet))
    except:
        print('distancetoatom: aborting - input error!')
        return False

    # round origin
    ori = [round(x,decimals) for x in ori]

    # report?
    if filename in ['', '0', 'False', 'None']:
        filename=False
    else:
        try:
            report=open(filename,'w') # file for writing
        except:
            print('distancetoatom: Unable to open report file! - Aborting!')
            return False

    # temporary name for pseudoatom
    tempname = cmd.get_unused_name('temp_name')
    tempsel = cmd.get_unused_name('temp_sel')

    #origin
    cmd.pseudoatom(object=tempname, resi=1, pos=ori)

    # select atoms within cutoff
    cmd.select(tempsel, '(%s around %f) and (%s) and state %d' %(tempname, cutoff, selection, state))
    cmd.delete(tempname)

    # single atom ori and part of selection
    # avoid double reporting!
    single_atom_ori=False
    try:
        if cmd.count_atoms('(%s) and (%s) and (%s)'%(selection, tempsel, origin))==1:
            single_atom_ori=True
    except: pass
    # pass= coordinates or multi-atom or single, not selected --> report ori

    # atom list
    stored.temp=[]
    cmd.iterate(tempsel, 'stored.temp.append([model, segi, chain, resn, resi, name, alt])')

    # custom properties? # conditional disabling
    if (property_name==''): property_name=False
    if ((cmd.get_version()[1]<1.7) and (property_name not in ['b','q'])):
        property_name=False

    # calculate the distances, creating list
    distance_list=[]
    if (not single_atom_ori):
        distance_list.append(['ORIGIN: '+str(origin), ori[0], ori[1], ori[2], 0.0])
    for atom in stored.temp:
        atom_name = ('/%s/%s/%s/%s`%s/%s`%s'%(atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6]))
        atom_xyz = [round(x, decimals) for x in cmd.get_atom_coords(atom_name)]
        atom_dist = round(cpv.distance(ori, atom_xyz), decimals)
        distance_list.append([atom_name,atom_xyz[0],atom_xyz[1],atom_xyz[2], atom_dist])
        # create property with distance (can be used for coloring, labeling etc)
        if property_name:
            try:
                cmd.alter(atom_name, '%s=%f'%(property_name, atom_dist))
            except:
                # I'm not sure alter raises exceptions if the property is illegal
                property_name=False

    # sort list, if selected
    if sort>0: distance_list.sort(key=lambda dist: dist[4])
    elif sort<0: distance_list.sort(key=lambda dist: dist[4], reverse=True)
    # add header
    distance_list=[['Atom Macro ID',
                    'x-coord',
                    'y-coord',
                    'z-coord',
                    'distance_to_origin']
                ]+distance_list

    if ((not quiet) and (filename)):
        # Hijack stdout to print to file and console
        class logboth(object):
            def __init__(self, *files):
                self.files = files
            def write(self, obj):
                for f in self.files:
                    f.write(obj)
        originalstdout = sys.stdout
        sys.stdout = logboth(sys.stdout, report)

    for entry in distance_list:
        if coordinates:
            output= '%s, %s, %s, %s, %s' %(entry[0],entry[1],entry[2],entry[3],entry[4]) #csv style
        else:
            output= '%s, %s' %(entry[0],entry[4]) #csv style
        if (not quiet):
            print(output)
        elif filename:
            report.write(output+'\n')

    # restore regular stdout
    if ((not quiet) and (filename)): sys.stdout = originalstdout
    # close file
    if filename: report.close()

    if (not quiet):
        if property_name: print('Distances saved to property: %s' %str(property_name))
        else: print('Distances NOT saved to property (illegal custom property)')

    # remove temp. selection
    cmd.delete(tempsel)

    # return list for potential use:
    if coordinates:
        if len(distance_list)>2: # prevents crash if list is otherwise empty
            distance_list2=list(map(distance_list.__getitem__, [1,4]))
            return distance_list2
        else: return distance_list
    else:
        return distance_list
################################################################################
cmd.extend( 'distancetoatom', distancetoatom );
cmd.auto_arg[0]['distancetoatom']=[lambda: cmd.Shortcut(['all', '[0,0,0]', 'pk1']), 'origin=', ',']
cmd.auto_arg[1]['distancetoatom']=[lambda: cmd.Shortcut(['7']), 'cutoff=', ',']
cmd.auto_arg[2]['distancetoatom']=[lambda: cmd.Shortcut(['None','distancetoatom_report.txt']), 'filename=', ',']
cmd.auto_arg[3]['distancetoatom']=[lambda: cmd.Shortcut(['selection=all, state=0, property_name=p.dist, coordinates=0, decimals=3, sort=1']), 'remaining (defaults)...', '']
################################################################################
