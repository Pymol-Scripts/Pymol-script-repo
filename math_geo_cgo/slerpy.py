################################################################################
#slerpy.py - Command definition routines for slerpy
#Copyright (C) 2006 Joel Bard
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#################################################################################
 
from pymol import cmd
import movs
 
def readViews( filename ):
#obsolete
    vfile = open( filename, 'r')
    views = [] #a list of views each of which is a list of 18 numbers
    vstrings = vfile.readlines()
    for vstring in vstrings:
        vals = vstring.split()
        view = [float(v) for v in vals]
        views.append( view )
    vfile.close()
    return views
 
def readFrames( filename ):
#obsolete
    ffile = open( filename, 'r' )
    frames = []
    fstrings = ffile.readlines()
    for fstring in fstrings:
        frames.append( int(fstring) )
    ffile.close()
    return frames
 
def readActions( filename ):
#obsolete
    #actions are stored in file where
    #for each line, the first 4 chars are the view index
    #associated with the action and the rest of the
    #line is the pymol command to be executed
    #upon reading, a dictionary is returned
    #with view indices as keys and actions as values
    actions = {}
    try:
        afile = open( filename, 'r' )
    except:
        print "No actions for this project"
        return actions
    astrings = afile.readlines()
    for astring in astrings:
        try:
            aindex = int(astring[:4])
            action = astring[4:]
            actions[ aindex ] = action[:-1]
        except:
            print "empty line"
    afile.close()
    return actions
 
def readModels( filename ):
#obsolete
    models = {}
    try:
        mfile = open( filename, 'r' )
    except:
        print "No models for this project"
        return models
    mstrings = mfile.readlines()
    for mstring in mstrings:
        try:
            mindex = int(mstring[:4])
            model = mstring[4:]
            models[ mindex ] = model[:-1]
        except:
            print "empty line"
    mfile.close()
    return models
 
def readSettings( filename ):
#obsolete
    settings = {}
    try:
        sfile = open( filename, 'r' )
    except:
        print "No settings for this project"
        return settings
    sstrings = sfile.readlines()
    for sstring in sstrings:
        try:
            sindex = int(sstring[:4])
            scommas = sstring[4:]
            settingName,selection,startVal,endVal = scommas.split(',')
            setting = [settingName,selection,float(startVal),float(endVal)]
            settings[sindex] = setting
        except:
            print "unable to parse setting"
    sfile.close()
    return settings
 
def readScenes( filename ):
#obsolete
    global scene_counter
 
    scene_counter = 0
    scenes = {}
    try:
        sfile = open( filename, 'r' )
    except:
        print "No scenes file for this project"
        return scenes
    sstrings = sfile.readlines()
    for sstring in sstrings:
        try:
            sindex = int(sstring[:4])
            scene = sstring[4:]
            scenes[ sindex ] = scene[:-1]
            scene_counter += 1
            #print "reading scene", sstring, sindex, scene
        except:
            print "empty line"
    sfile.close()
    return scenes
 
def read_all( fileroot ):
#obsolete in favor of readKeyViewFile
    global views
    global frames
    global actions
    global cur_view
    global cur_index
    global scenes
    global models
    global settings
 
    views = readViews( fileroot+".txt" )
    frames = readFrames( fileroot+".frm")
    actions = readActions( fileroot+".act")
    scenes = readScenes( fileroot+".scn")
    models = readModels( fileroot+".mod")
    settings = readSettings( fileroot+".set")
    cur_view = views[0]
    cur_index = 0
    show_cur()
 
def print_view( view ):
    for i in range(0,6):
        for j in range(0,3):
            print "%12.6f"% view[ 3*i+j ] ,
        print
 
def show_next():
    global cur_index
    global cur_view
    global views
    if( cur_index == len( views )-1 ):
        print "No more views."
        return
    cur_index += 1
    cur_view = views[ cur_index ]
    cmd.set_view(cur_view)
    print "Matrix: "
    print_view( cur_view )
    print "Showing view number ",cur_index
 
def show_prev():
    global cur_index
    global cur_view
    global views
    if( cur_index == 0 ):
        print "No more views."
        return
    cur_index -= 1
    cur_view = views[ cur_index ]
    cmd.set_view(cur_view)
    print "Matrix: "
    print_view( cur_view )
    print "Showing view number ",cur_index
 
def show_cur():
    global cur_index
    global cur_view
    global views
    cur_view = views[ cur_index ]
    cmd.set_view(cur_view)
    print "Matrix: "
    print_view( cur_view )
    print "Showing view number ",cur_index
 
def go_to_view( arg=0 ):
    global cur_index
    global cur_view
    global views
    n = int( arg )
    if n < 0 or n > len(views)-1:
        print "Index out of range."
        return
    cur_index = n
    cur_view = views[n]
    cmd.set_view( cur_view )
    print "Matrix: "
    print_view( cur_view )
    print "Showing view number ", cur_index
 
def insert_current():
    #insert the current view into the list after the view
    #in views[cur_index]
    #set frames to default
    global cur_index
    global cur_view
    global views
    global frames
    global actions
    global scenes
    global settings
    global models
    global fades
 
    cur_index += 1
    cur_view = cmd.get_view()
    views.insert( cur_index, [cv for cv in cur_view] )
    frames.insert( cur_index, 50 )
 
    #deal with actions dictionary
    actions = incKeyAbove( actions, cur_index )
    scenes = incKeyAbove( scenes, cur_index )
    settings = incKeyAbove( settings, cur_index )
    models = incKeyAbove( models, cur_index )
    fades = incKeyAbove( fades, cur_index )
 
    print "New view:"
    print_view( cur_view )
    print "Inserted view with index", cur_index, "and a 50 frame transition"
 
def append_view():
    global views
    global frames
    global cur_index
    global cur_view
 
    cur_index = len(views)
    cur_view = cmd.get_view()
    views.append( [cv for cv in cur_view] )
    frames.append( 50 )
 
    print "New view: "
    print_view( cur_view )
    print "Appended view with index", cur_index, "and a 50 frame transition"
    print "The current view is", cur_index
 
def incKeyAbove( dict, index ):
    tempDict = {}
    for key, val in dict.iteritems():
        if key >= index:
            newkey = key + 1
            tempDict[newkey] = val
        else:
            tempDict[key] = val
    return tempDict
 
def decKeyAbove( dict, index ):
    tempDict = {}
    for key, val in dict.iteritems():
        if key > index:
            newkey = key - 1
            tempDict[newkey] = val
        else:
            tempDict[key] = val
    return tempDict
 
def delete_current():
    #remove the current view from the list
    #show the previous view
    global cur_index
    global cur_view
    global views
    global actions
    global scenes
    global settings
    global models
    global frames
    global fades
 
    del views[cur_index]
    del frames[cur_index]
    if actions.has_key( cur_index ):
        del actions[cur_index]
    if scenes.has_key( cur_index ):
        del scenes[cur_index]
    if settings.has_key( cur_index ):
        del settings[cur_index]
 
    #deal with dictionaries
    actions = decKeyAbove( actions, cur_index )
    scenes = decKeyAbove( scenes, cur_index )
    settings = decKeyAbove( settings, cur_index )
    models = decKeyAbove( models, cur_index )                               
    fades = decKeyAbove( fades, cur_index )
 
    print "View number",cur_index,"deleted."
    if cur_index > 0:
        cur_index -= 1
    cur_view = views[cur_index]
    cmd.set_view( cur_view )
    print "Current view is number",cur_index
 
def delete_settings():
    global settings
    global cur_index
    del settings[cur_index]
 
def replace_current():
    global cur_index
    global cur_view
    global views
    cur_view = cmd.get_view()
    views[cur_index] = [cv for cv in cur_view]
 
def insert_scene():
    global views
    global actions
    global settings
    global frames
    global cur_index
    global cur_view
    global scenes
    global scene_counter
    global models
    global fades
 
    cur_index += 1
    cur_view = cmd.get_view()
    views.insert( cur_index, [cv for cv in cur_view] )
    frames.insert( cur_index, 50 )
 
    #deal with dictionaries
    actions = incKeyAbove( actions, cur_index )
    scenes = incKeyAbove( scenes, cur_index )
    settings = incKeyAbove( settings, cur_index )
    models = incKeyAbove( models, cur_index )
    fades = incKeyAbove( fades, cur_index )
 
    #this stuff has to be done after the above
    #find a free scene name
    i = 1
    found = 1
    while found == 1:
        found = 0
        for sname in scenes.values():
            print "|"+sname+"|"
            if sname == "slerpy_"+str(i):
                found = 1
                break
        if found == 1:
            i += 1
        else:
            break
    newname = "slerpy_"+str(i)
 
    scene_counter += 1
    cmd.scene( newname, "store" )
    scenes[ cur_index ] = newname
    actions[ cur_index ] = "scene "+newname
 
    print "New view:"
    print_view( cur_view )
    print "Inserted view with index", cur_index, "and a 50 frame transition"
    print "Added scene",newname
 
def write_views( filename ):
#deprecated in favor of key files
    global views
    global frames
    global scenes
    global actions
    global settings
    global models
    global fades
 
    viewfile = open( filename+".txt", 'w')
    for view in views:
        for v in view:
            viewfile.write( str(v) + " " )
        viewfile.write('\n')
    viewfile.close()
 
    framefile = open( filename+".frm", 'w' )
    for frame in frames:
        framefile.write( str( frame )+'\n' )
    framefile.close()
 
    actionfile = open( filename+".act", 'w' )
    for key,action in actions.iteritems():
        keystring = str( key )
        actionfile.write( keystring.rjust(4)+action + '\n' )
    actionfile.close()
 
    scenefile = open( filename+".scn", 'w' )
    for key,scene in scenes.iteritems():
        keystring = str( key )
        scenefile.write( keystring.rjust(4)+scene + '\n' )
    scenefile.close()
 
    modelfile = open( filename+".mod", 'w' )
    for key,model in models.iteritems():
        keystring = str( key )
        modelfile.write( keystring.rjust(4)+model + '\n' )
    modelfile.close()
 
    settingsfile = open( filename+".set", 'w' )
    for key,setting in settings.iteritems():
        keystring = str( key )
        settingName, selection, startVal, endVal = setting
        settingsfile.write( "%s%s,%s,%f,%f\n" % (keystring.rjust(4), settingName, selection, startVal, endVal))
    settingsfile.close()
    cmd.save( filename+".pse" )
 
    print "Wrote files", filename+".txt,",filename+".frm,",filename+".pse, and",filename+".act"
 
def writeKeyViewFile( filename ):
    global views
    global frames
    global scenes
    global actions
    global settings
    global models
    global fades
 
    keyviewfile = open( filename + ".key", 'w' )
    i=0
    for view in views:
        keyviewfile.write( "VIEW: %4d " % i )
        for v in view:
            keyviewfile.write( str(v) + " " )
        keyviewfile.write('\n')
        keyviewfile.write( "FRAMES: %d\n" % frames[i] )
        if actions.has_key( i ):
            keyviewfile.write( "ACTIONS: %s\n" % actions[i] )
        if scenes.has_key( i ):
            keyviewfile.write( "SCENES: %s\n" % scenes[i] )
        if models.has_key( i ):
            keyviewfile.write( "MODELS: %s\n" % models[i] )
        if settings.has_key( i ):
            settingName, selection, startVal, endVal = settings[i]
            keyviewfile.write( "SETTINGS: %s, %s, %f, %f\n" % (settingName, selection, startVal, endVal))
        if fades.has_key( i ):
            startVisSelection, endVisSelection, sticksOnly = fades[i]
            keyviewfile.write( "FADES: %s, %s, %d\n" % (startVisSelection, endVisSelection, sticksOnly))            
        i += 1
        keyviewfile.write("\n")
    keyviewfile.close()
    cmd.save( filename + ".pse" )
    print "Wrote files " , filename + ".key", filename + ".pse"
 
def readKeyViewFile( filename ):
    global views
    global frames
    global scenes
    global actions
    global settings
    global models
    global fades
    global scene_counter
 
    views = []
    frames = []
    actions = {}
    scenes = {}
    models = {}
    settings = {}
    fades = {}
    scene_counter = 0
 
    if filename[-4:] == ".key": filename = filename[:-4]                                                                                                                  
    keyviewfile = open( filename + ".key", 'r' )
    viewstrings = keyviewfile.readlines()
    keyviewfile.close()
    viewindex = 0
    for line in viewstrings:
        if line.startswith("VIEW: "):
            viewindex = int( line[6:10] )
            vals = line[10:].split()
            view = [float(v) for v in vals]
            views.append( view )
        if line.startswith("FRAMES: "):
            frames.append( int( line[8:] ) )
        if line.startswith("ACTIONS: "):
            actions[ viewindex ] = line[9:-1]
        if line.startswith("SCENES: "):
            scenes[ viewindex ] = line[8:-1]
            scene_counter += 1
        if line.startswith("MODELS: "):
            models[ viewindex ] = line[8:-1]
        if line.startswith("SETTINGS: "):
            settingName,selection,startVal,endVal = line[10:-1].split(',')
            settings[ viewindex ] = [settingName,selection,float(startVal),float(endVal)]
        if line.startswith( "FADES: " ):
            startVisSelection, endVisSelection, sticksOnly = line[7:-1].split(',')
            fades[ viewindex ] = [startVisSelection, endVisSelection, int(sticksOnly) ]
    cur_view = views[0]
    cur_index = 0
    show_cur()
 
def set_frames_current( arg=50 ):
    global frames
    global cur_index
    frames[cur_index] = int(arg)
 
def list_frames():
    global frames
    global views
    global actions
    global models
    global settings
 
    i=0
    f=0
    for view in views[:-1]:
        if actions.has_key( i ):
            a = actions[i]
        else:
            a = ""
        if models.has_key( i ):
            m = "States: " + models[i]
        else:
            m = ""
        if settings.has_key( i ):
            settingName, selection, startVal, endVal = settings[i]
            s = "Settings: %s %s %f %f" % (settingName, selection, startVal, endVal)
        else:
            s = ""
        print "View",i,"to",i+1,"Frames ",f,"to",f+frames[i],a,m,s
        f += frames[i]
        i += 1
 
def add_action_current( cmd ):
    global cur_index
    global actions
    actions[cur_index] = cmd[1:-1] #strip off quotes
 
def append_action_current( cmd ):
    global cur_index
    global actions
    actions[cur_index] += ";" + cmd[1:-1]
 
def clear_action_current():
    global actions
    global cur_index
    del actions[cur_index]
 
def list_actions():
    global actions
    for i,a in actions.iteritems():
        print i,a
 
def morph_models( start_model, end_model ):
    global cur_index
    global frames
    global models
    models[cur_index] = "%s -%s" % (start_model, end_model)
    frames[cur_index] = abs(int(end_model) - int(start_model)) + 1
 
def interpolate_settings( setting, selection, startval, endval ):
    global cur_index
    global settings
    settingEntry = [setting, selection, float(startval), float(endval)]
    settings[cur_index] = settingEntry 
 
def crossfade( startVisSelection, endVisSelection, sticksOnly = 1 ):
#cross fade the specified objects, vary stick transparency only if stickOnly=1
    global cur_index
    global fades
    fades[cur_index] = [startVisSelection, endVisSelection, int(sticksOnly) ]
 
def setup_view( index ):
    for i in range( 0, int(index)+1 ):
        if actions.has_key(i):
            print "Executing %s from actions %d" % (actions[i],i)
            cmd.do( actions[i] )
        if settings.has_key(i):
            settingName, selection, startVal, endVal = settings[i]
            action = "set %s, %f, %s;" % (settingName, endVal, selection)
            print "Executing %s from settings %d" % (action,i)
            cmd.do( action )
        if fades.has_key(i):
            startVisSelection, endVisSelection, sticksOnly = fades[i]
            action = "set stick_transparency, 0, %s; set stick_transparency, 1, %s;" % (endVisSelection, startVisSelection)
            print "Executing %s from fades %d" % (action, i)
            cmd.do( action )
 
def show_transition(start_index=0, end_index=0):
    #show the transition from the current view to the next view
    global frames
    global views
    global cur_index
    global actions
    global models
    if( start_index == 0 and end_index == 0 ):
        if cur_index >= len(views)-1:
            print "Current view is last in sequence."
            return
        start_index=cur_index
        end_index=cur_index+1
    else:
        start_index = int(start_index)
        end_index = int(end_index)
    ftot = 0
    setcommand = ""
    i = start_index
    for nframes in frames[start_index:end_index]:
        #ftot += nframes
        if models.has_key(i):
            setcommand += " " + models[i] + " "
        else:
            setcommand += " 1 x%i" % nframes
        i += 1
 
#    cmd.mset("1 x%i" % ftot)
    cmd.mset( setcommand )
    start_frame = 1
    #first do all actions that happen up to this point to make sure
    #things look the way they should
    first_action = ""
    for i in range( 0, start_index ):
        if actions.has_key(i):
            first_action += actions[i] + ';'
            #print "Executing %s from actions %d" % (actions[i],i)
            #cmd.do( actions[i] )
        if settings.has_key(i):
            settingName, selection, startVal, endVal = settings[i]
            action = "set %s, %f, %s;" % (settingName, endVal, selection)
            first_action += action
            #print "Executing %s from settings %d" % (action,i)
            #cmd.do( action )
        if fades.has_key(i):
            startVisSelection, endVisSelection, sticksOnly = fades[i]
            action = "set stick_transparency, 0, %s; set stick_transparency, 1, %s;" % (endVisSelection, startVisSelection)
            first_action += action
            #print "Executing %s from fades %d" % (action, i)
            #cmd.do( action )
    for i in range( start_index, end_index ):
        if settings.has_key(i):
            movs.animate_transition( views[i], views[i+1], frames[i], start_frame, settings[i] )
        elif fades.has_key(i):
            movs.animate_transition( views[i], views[i+1], frames[i], start_frame, fades[i] )
        else:
            movs.animate_transition( views[i], views[i+1], frames[i], start_frame )
        #add an action
        if start_frame == 1:
            mdo_cmd = first_action
            if actions.has_key(i):
                mdo_cmd += actions[i]+";"
            #mdo_cmd += "set_view("+str(views[i])+")"
            print mdo_cmd
            cmd.mdo(start_frame, mdo_cmd)
        elif actions.has_key(i):
            mdo_cmd = actions[i]+";set_view("+str(views[i])+")"
            cmd.mdo(start_frame, mdo_cmd)
            #print mdo_cmd
        start_frame += frames[i]
    cmd.frame(1)
    cmd.mplay()
 
def make_all():
    #make the whole movie
    global views
    global frames
    global models
 
    #first get total number of frames
    ftot = 0
    i = 0
    setcommand = ""
    for nframes in frames[:-1]:
        ftot += nframes
        if models.has_key(i):
            setcommand += " " + models[i] + " "
        else:
            setcommand += " 1 x%i" % nframes
        i += 1
 
    #initialize movie
    #cmd.mset("1 x%i" % ftot)
    #cmd.mset("1 x50 1 -30 30 x20")
    cmd.mset( setcommand )
 
    #loop through views
    start_view = views[0][:]
    i = 0
    first_frame = 1
    for view in views[1:]:
        end_view = view[:]
        if settings.has_key(i):
            movs.animate_transition( start_view, end_view, frames[i], first_frame, settings[i] )
        elif fades.has_key(i):
            movs.animate_transition( start_view, end_view, frames[i], first_frame, fades[i] )
        else:
            movs.animate_transition( start_view, end_view, frames[i], first_frame )
        #add an action
        if actions.has_key(i):
            mdo_cmd = actions[i]#+";set_view ("+str( views[i] )+")"
            print mdo_cmd
            cmd.mdo(first_frame, mdo_cmd)
        first_frame += frames[i]
        i += 1
        start_view = end_view[:]
    cmd.frame(1)
 
## views = readViews( "viewfile.txt" )
## frames = readFrames( "viewfile.frm" )
## actions = readActions( "viewfile.act" )
##print "Length ",len(views)
#for v in views:
 #   print v
#cur_view = views[0]
views = []
frames = []
models = {}
actions = {}
scenes = {}
settings = {}
fades = {}
scene_counter = 0
cur_index = -1
cmd.set( "scene_animation_duration","0" )
#cmd.set_view( cur_view )
 
cmd.extend("sn", show_next )
cmd.extend("sp", show_prev )
cmd.extend("sc", show_cur )
cmd.extend("sinsert", insert_current )
cmd.extend("sdelete", delete_current )
cmd.extend("sreplace", replace_current )
cmd.extend("sappend", append_view )
cmd.extend("sscene", insert_scene )
cmd.extend("sgo", go_to_view )
cmd.extend("sreadold", read_all )
cmd.extend("swriteold", write_views )
cmd.extend("slist", list_frames )
cmd.extend("ssetf", set_frames_current )
cmd.extend("sshow", show_transition )
cmd.extend("srecord", make_all )
cmd.extend("saction", add_action_current )
cmd.extend("sdelaction", clear_action_current )
cmd.extend("sdumpactions", list_actions )
cmd.extend("sappendaction", append_action_current )
cmd.extend("smorph", morph_models )
cmd.extend("sinterpsetting", interpolate_settings )
cmd.extend("sdeletesetting", delete_settings )
cmd.extend("scrossfade", crossfade )
cmd.extend("swrite", writeKeyViewFile )
cmd.extend("sread", readKeyViewFile )
cmd.extend("ssetupview", setup_view )
