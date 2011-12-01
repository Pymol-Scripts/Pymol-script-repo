from pymol import cmd
def testme():
    print("Hello world")
cmd.extend("testme",testme)
