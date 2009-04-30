#
# apropos.py 
# Author: Ezequiel Panepucci
# Date: 2006-07-20
#
from pymol import cmd
import re
 
def apropos(regexp=''):
    '''
DESCRIPTION
        "apropos" searches through the documentation of all currently 
        defined commands and lists those commands for which the keyword
        is either contained in the documentation or matches the command
        name itself.
 
        If an appropriate "DESCRIPTION" section is provided in the documentation
        of the command, the first 80 characters are listed as a summary.
 
USAGE
        apropos [keyword or regexp]
 
EXAMPLE
        apropos fit
 
###EXACT MATCH FOR: fit ==> try 'help fit' at the prompt.
 
###The following commands are NOT documented.
 
      vdw_fit
 
###The following commands are documented.  'help command' 
 
          fit : "fit" superimposes the model in the first selection on to the model
    intra_fit : "intra_fit" fits all states of an object to an atom selection
          rms : "rms" computes a RMS fit between two atom selections, but does not
     pair_fit : "pair_fit" fits a set of atom pairs between two models.  Each atom
intra_rms_cur : "intra_rms_cur" calculates rms values for all states of an object
     commands : >>>>>>>> Ooopsie, no DESCRIPTION found for this command!!! <<<<<<
         zoom : "zoom" scales and translates the window and the origin to cover the
    intra_rms : "intra_rms" calculates rms fit values for all states of an object
        align : "align" performs a sequence alignment followed by a structural
      rms_cur : "rms_cur" computes the RMS difference between two atom
      fitting : "fitting" allows the superpositioning of object1 onto object2 using
 
SEE ALSO
    grepset(www.pymolwiki.org), Python re module
    '''
    cmd.set("text","1",quiet=1)
 
    count=0
    docre = re.compile(regexp, re.MULTILINE | re.IGNORECASE)
    cmdre = re.compile(regexp, re.IGNORECASE)
 
    matches_with_help = []
    matches_without_help = []
 
    maxcclen=0
    for cc in cmd.keyword:
        if cc == regexp:
            print '\n###EXACT MATCH FOR: %s ==> try \'help %s\' at the prompt.' % (cc,cc)
 
        doc = cmd.keyword[cc][0].__doc__
 
        if doc == None:
            if re.search(regexp, cc, re.IGNORECASE):
                count += 1
                matches_without_help.append(cc)
            continue
 
        if re.search(regexp, doc, re.MULTILINE | re.IGNORECASE):
            count += 1
            if len(cc) > maxcclen:
                maxcclen = len(cc)
 
            docmatches = re.match(r"""^\s+DESCRIPTION\s+(.{0,80})\S*""", doc, re.IGNORECASE)
            if docmatches == None:
                desc = '>>>>>>>> Ooopsie, no DESCRIPTION found for this command!!! <<<<<<'
            else:
                desc = docmatches.groups()[0]
            matches_with_help.append( (cc, desc ) )
 
 
    if len(matches_without_help) > 0:
        fmt = '%' + str(maxcclen) + 's' # get the width of the lenghtiest command in a format
        print '\n###The following commands are NOT documented.\n'
        for cc in matches_without_help:
            print fmt % (cc)
 
    if len(matches_with_help) > 0:
        fmt = '%' + str(maxcclen) + 's : %s' # get the width of the lenghtiest command in a format
        print '\n###The following commands are documented.  \'help command\' \n'
        for (cc,desc) in matches_with_help:
            print fmt % (cc,desc)
 
cmd.extend('apropos',apropos)
