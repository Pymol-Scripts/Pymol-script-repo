#########################################################################
#
# Date: May 2003 Authors: Sophie COON
# 

#########################################################################
# This is a fix to use the askdirectory widget developped by written by
# Fredrik Lundh but only available for Python2.2 and higher.
#########################################################################

from tkCommonDialog import Dialog

class Directory(Dialog):
    "Ask for a directory"

    command = "tk_chooseDirectory"

    def _fixresult(self, widget, result):
        self.widget = widget
        if result:
            # keep directory until next time
            self.options["initialdir"] = result
        self.directory = result # compatibility
        return result

class CreateDirectory(Dialog):
    command = "tk_getSaveFile"
    def _fixresult(self, widget, result):
        if result:
            # keep directory until next time
            self.options["initialdir"] = result
        self.directory = result # compatibility
        return result


def askdirectory (**options):
    "Ask for a directory, and return the file name"
    return apply(Directory, (), options).show()

def createdirectory(**options):
    "Ask for a directory, and return the file name"
    return apply(CreateDirectory, (), options).show()
    
