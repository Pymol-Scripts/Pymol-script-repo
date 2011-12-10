# $Header: /opt/cvs/python/packages/share1.5/mglutil/splashregister/license.py,v 1.15 2010/08/19 20:14:31 sargis Exp $
# $Id: license.py,v 1.15 2010/08/19 20:14:31 sargis Exp $
#

import Tkinter

tk_root = Tkinter.Tk()
tk_root.title("Commercial Usage")
txt = """
 The software component for computing molecular surfaces (MSMS) 
 is not free for commercial usage. If you plan to use MSMS for commercial 
 research please contact sanner@scripps.edu

 Some software components such at the volume rendering and 
 isocontouring were developed at UT Austin.

 If you publish scientific results generated using this software 
 please cite the appropriate software components.
 A list of papers is provided under Help -> Citation Information 
 menu in PMV and ADT. 
"""
Tkinter.Label(tk_root, text=txt, justify=Tkinter.LEFT).pack()
Tkinter.Button(tk_root, text="OK", command=tk_root.quit).pack()
tk_root.mainloop()
