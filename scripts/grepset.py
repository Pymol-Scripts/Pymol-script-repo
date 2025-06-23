from pymol import cmd
import re
import pymol.setting


def grepset(regexp=''):
    '''
DESCRIPTION
    "grepset" greps through the list of settings using a python
    regular expression as defined in the 're' module.
    It returns a list of settings/values matching the regexp.
    No regexp returns every setting.

USAGE
    grepset [regexp]

EXAMPLE
    grepset line
    grepset ray
    grepset (^line|color$)

SEE ALSO
        Python re module
    '''

    count = 0
    regexp = re.compile(regexp)
    matches = []
    for a in pymol.setting.get_index_list():
        setting = pymol.setting._get_name(a)
        if regexp.search(setting):
            count += 1
            matches.append((setting, cmd.get_setting_text(a, '', -1)))
    # max length of the setting names that matched
    maxlen = max([len(s[0]) for s in matches] + [0])
    fmt = "%%-%ds : %%s" % (maxlen,)
    for setting in matches:
        print((fmt % setting))
    print(('%d settings matched' % (count,)))
    cmd.set('text', 1)
cmd.extend('grepset', grepset)
