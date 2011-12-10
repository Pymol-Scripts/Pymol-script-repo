#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/moleculeParser.py,v 1.9 2008/09/11 17:00:57 sargis Exp $
#
# $Id: moleculeParser.py,v 1.9 2008/09/11 17:00:57 sargis Exp $
#

from MolKit.molecule import Atom
import warnings, os
from mglutil.util.misc import ensureFontCase

class MoleculeParser:
    def __init__(self, filename=None, allLines=None):
        """Supply the filename for reading by readFile, or
        supply the lines directly via allLines
        """
        self.filename = filename
        self.allLines = allLines #stores all lines from file

        
    def readFile(self):
        f = open(self.filename)
        self.allLines = f.readlines()
        if len(self.allLines)==1:
            # this file probably has \r instead or \n
            self.allLines = self.allLines[0].split('\r')
            warnings.warn('Only 1 line read from PDB file, splitting on \r')
        f.close()
        import string
        self.allLines = filter( lambda x,s=string.strip: len(s(x)),
                                self.allLines )

    
    def viewSource(self):
        import Tkinter, Pmw
        root = Tkinter.Toplevel()
        root.title(self.filename)
        self.st = Pmw.ScrolledText(root)
        self.st.pack(fill = 'both', expand=1)
        
        self.st._textbox.configure(bg='white', font=(ensureFontCase('Courier'), '10'))
        txt = ''
        for line in self.allLines:
            txt += ''.join(line)
        self.st.setvalue(txt)