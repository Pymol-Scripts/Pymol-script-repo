from pymol import cmd
def testme():
    print("Hello world")
    print("Hello changes")
cmd.extend("testme",testme)
