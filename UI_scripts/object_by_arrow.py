from pymol import cmd
 
def move_down():
        enabled_objs = cmd.get_names("object",enabled_only=1)
        all_objs = cmd.get_names("object",enabled_only=0)
        for obj in enabled_objs:
                cmd.disable(obj)
                last_obj=obj
                for i in range(0,len(all_objs)):
                        if all_objs[i] == obj:
                                if i+1 >= len(all_objs):
                                        cmd.enable( all_objs[0] )
                                else:
                                        cmd.enable( all_objs[i+1] )
        cmd.orient
def move_up():
        enabled_objs = cmd.get_names("object",enabled_only=1)
        all_objs = cmd.get_names("object",enabled_only=0)
        for obj in enabled_objs:
                cmd.disable(obj)
                last_obj=obj
                for i in range(0,len(all_objs)):
                        if all_objs[i] == obj:
                                if i-1 < 0:
                                        cmd.enable( all_objs[-1] )
                                else:
                                        cmd.enable( all_objs[i-1] )
        cmd.orient
 
cmd.set_key('left', move_up)
cmd.set_key('right', move_down)
