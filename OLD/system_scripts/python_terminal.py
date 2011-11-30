from pymol import cmd
 
def parse(command):
	exec command
 
cmd.extend('py', parse)

