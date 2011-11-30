## This is just a quick hack. For something more meaty see;
## http://arcib.dowling.edu/sbevsl/
 
## Version 0.0.00-000/1
 
 
## Turn off the virtual_trackball
cmd.set("virtual_trackball", "off")
 
 
## spacefill
def spacefill(p1=''):
    if(p1=='off'):
        cmd.hide("spheres")
    elif(p1==''):
        cmd.show("spheres")
    else:
        print("feh!")
cmd.extend("spacefill", spacefill)
 
## cartoon
def cartoon(p1=''):
    if(p1=='off'):
        cmd.hide("cartoon")
    elif(p1==''):
        cmd.show("cartoon")
    else:
        print("feh!")
cmd.extend("cartoon", cartoon)
 
## wireframe
def wireframe(p1=''):
    if(p1=='off'):
        cmd.hide("lines")
    elif(p1==''):
        cmd.show("lines")
    else:
        print("feh!")
cmd.extend("wireframe", wireframe)
 
 
## exit
def exit():
    cmd.quit()
cmd.extend("exit", exit)

