from pymol import cmd
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
 
SEE ALSO
   Python re module
   '''
 
   from re import compile
 
   count=0
   regexp=compile(regexp)
   for a in pymol.setting.get_index_list():
      setting=pymol.setting._get_name(a)
      if regexp.search(setting):
         count = count + 1
         print '%-30s %s' % (setting, cmd.get_setting_text(a,'',-1))
 
   print '%d settings matched' % count
cmd.extend('grepset',grepset)
