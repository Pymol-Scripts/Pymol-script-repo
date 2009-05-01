#script contributed by Philippe Garteiser; garteiserp@omrf.ouhsc.edu
from pymol import cmd
 
def resicolor(selection='all'):
 
    '''USAGE: resicolor <selection>
    colors all or the given selection with arbitrary
    coloring scheme.
    '''
    cmd.select ('calcium','resn ca or resn cal')
    cmd.select ('acid','resn asp or resn glu or resn cgu')
    cmd.select ('basic','resn arg or resn lys or resn his')
    cmd.select ('nonpolar','resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
    cmd.select ('polar','resn ser or resn thr or resn asn or resn gln or resn tyr')
    cmd.select ('cys','resn cys or resn cyx')
    cmd.select ('backbone','name ca or name n or name c or name o')
    cmd.select ('none')
 
    print selection
    code={'acid'    :  'red'    ,
          'basic'   :  'blue'   ,
          'nonpolar':  'orange' ,
          'polar'   :  'green'  ,
          'cys'     :  'yellow'}
    cmd.select ('none')
    for elem in code:
        line='color '+code[elem]+','+elem+'&'+selection
        print line
        cmd.do (line)
    word='color white,backbone &'+selection
    print word
    cmd.do (word)                  #Used to be in code, but looks like
                                   #dictionnaries are accessed at random
    cmd.hide ('everything','resn HOH')
 
cmd.extend ('resicolor',resicolor)

