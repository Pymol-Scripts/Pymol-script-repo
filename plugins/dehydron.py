'''
Dehydron: A dehydron calculator plugin for PyMOL
Version: 1.5
Described at PyMOL wiki:
http://www.pymolwiki.org/index.php/dehydron

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date    : January 2012
License: GNU General Public License
Acknowledgement: The H-bond detection code is based on the list_mc_hbonds.py
script from Robert L. Campbell http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
'''

import Tkinter
from Tkinter import *
import Pmw
from pymol import cmd
from pymol import stored


def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'dehydron',
                            label = 'Dehydron',
                            command = lambda : mainDialog())


def mainDialog():
    """ Creates the GUI"""

    def get_dehydrons():
        angle_range = float(angle_value.get())
        max_distance = float(dist_cutoff_value.get())
        desolv = float(desolv_sphere.get())
        min_wrappers = float(min_value.get())
        selection = sel_value.get()
        dehydron(selection, angle_range, max_distance, desolv, min_wrappers)

    master = Tkinter.Tk()
    master.title(' dehydron ')
    w = Tkinter.Label(master, text = 'dehydron calculator\nOsvaldo Martin - aloctavodia@gmail.com',
                                background = '#000000',
                                foreground = '#cecece',
                                #pady = 20,
                                )
    w.pack(expand=1, fill = 'both', padx = 4, pady = 4)

    Pmw.initialise(master)
    nb = Pmw.NoteBook(master, hull_width = 420, hull_height=280)
    p1 = nb.add('Main')
    p2 = nb.add('About')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ Main TAB #################################
### hydrogen bond settings
    group = Pmw.Group(p1,tag_text='Hydrogen bond Settings')
    group.pack(fill='x', expand=1, padx=20, pady=1)
    Label(group.interior(), text='angle range').grid(row=2, column=0)
    angle_value = StringVar(master=group.interior())
    angle_value.set(40)
    entry_angle = Entry(group.interior(),textvariable=angle_value, width=10)
    entry_angle.grid(row=2, column=1)
    entry_angle.configure(state='normal')
    entry_angle.update()
    Label(group.interior(), text='max distance').grid(row=3, column=0)
    dist_cutoff_value = StringVar(master=group.interior())
    dist_cutoff_value.set(3.5)
    entry_dist = Entry(group.interior(),textvariable=dist_cutoff_value, width=10)
    entry_dist.grid(row=3, column=1)
    entry_dist.configure(state='normal')
    entry_dist.update()
### dehydron settings
    group = Pmw.Group(p1,tag_text='Dehydron Settings')
    group.pack(fill='x', expand=1, padx=20, pady=5)
    Label(group.interior(), text='desolvatation sphere radius').grid(row=2, column=2)
    desolv_sphere = StringVar(master=group.interior())
    desolv_sphere.set(6.5)
    entry_desolv=Entry(group.interior(),textvariable=desolv_sphere, width=10)
    entry_desolv.grid(row=2, column=3)
    entry_desolv.configure(state='normal')
    entry_desolv.update()
    Label(group.interior(), text='minimum wrappers').grid(row=3, column=2)
    min_value = StringVar(master=group.interior())
    min_value.set(19)
    entry_min_value=Entry(group.interior(),textvariable=min_value, width=10)
    entry_min_value.grid(row=3, column=3)
    entry_min_value.configure(state='normal')
    entry_min_value.update()
### selection settings
    group = Pmw.Group(p1,tag_text='Selection')
    group.pack(fill='x', expand=1, padx=20, pady=5)
    Label(group.interior(), text='object selection').grid(row=4, column=2)
    sel_value = StringVar(master=group.interior())
    sel_value.set('all')
    entry_sel_value=Entry(group.interior(),textvariable=sel_value, width=10)
    entry_sel_value.grid(row=4, column=3)
    entry_sel_value.configure(state='normal')
    entry_sel_value.update()
### submit
    Button(p1, text="Calculate", command=get_dehydrons).pack(side=BOTTOM)
############################ About TAB #################################
    group = Pmw.Group(p2, tag_text='About dehydron plug-in')
    group.pack(fill = 'both', expand=1, padx = 5, pady = 5)
    text =u"""For a brief introduction to the dehydron concept, you could
read http://en.wikipedia.org/wiki/dehydron

Citation for this plugin:
Martin O.A.; dehydron calculator plugin for PyMOL,
2012. IMASL-CONICET.

Citation for PyMOL may be found here:
http://pymol.sourceforge.net/faq.html#CITE

Citation for dehydrons (I think these could be used):
Fern\u00E1ndez A. and Scott R.; "Adherence of packing
defects in soluble proteins", Phys. Rev. Lett. 91, 018102
(2003).

Fern\u00E1ndez A., Rogale K., Scott R. and Scheraga H.A.;
"Inhibitor design by wrapping packing defects in HIV-1
proteins", PNAS, 101, 11640-45 (2004).

Fern\u00E1ndez A. "Transformative Concepts for Drug Design:
Target Wrapping" (ISBN 978-3642117916),
Springer-Verlag, Berlin, Heidelberg (2010).
"""
    #
    # Add this as text in a scrollable panel.
    # Code based on Caver plugin
    # http://loschmidt.chemi.muni.cz/caver/index.php
    #
    interior_frame = Frame(group.interior())
    bar = Scrollbar(interior_frame)
    text_holder = Text(interior_frame, yscrollcommand=bar.set, foreground="#cecece",background="#000000",font="Times 12")
    bar.config(command=text_holder.yview)
    text_holder.insert(END,text)
    text_holder.pack(side=LEFT,expand="yes",fill="both")
    bar.pack(side=LEFT,expand="yes",fill="y")
    interior_frame.pack(expand="yes",fill="both")

    master.mainloop()


def dehydron(selection='all', angle_range=40, max_distance=3.5, desolv=6.5, min_wrappers=19, quiet=0):
    '''
DESCRIPTION

    dehydron calculator

USAGE

    dehydron [ selection [, angle_range [, max_distance [, desolv [, min_wrappers ]]]]]
    '''

    angle, max_distance = float(angle_range), float(max_distance)
    desolv, min_wrappers = float(desolv), int(min_wrappers)
    quiet = int(quiet)

    name = 'DH_%s' % selection
    cmd.delete(name)

    selection = '((%s) and polymer)' % (selection)
    hb = cmd.find_pairs("((byres "+selection+") and n. n)","((byres "+selection+") and n. o)",mode=1,cutoff=max_distance,angle=angle_range)

    if not quiet:
        hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))
        print "------------------------------------------------"
        print "----------------Dehydron Results----------------"
        print "------------------------------------------------"
        print "    Donor      |    Aceptor    |  "
        print "Chain Residue  | Chain Residue | # wrappers"

    cmd.select('_nonpolar', 'not (solvent or hydro or (elem N+O) extend 1)', 0)

    sel = []
    for pairs in hb:
        wrappers = cmd.count_atoms('((%s and _nonpolar) within %f of byca (%s`%d %s`%d))' % 
                ((pairs[0][0], desolv) + pairs[0] + pairs[1]))
        if wrappers < min_wrappers:
            cmd.distance(name, pairs[0], pairs[1])
            if not quiet:
                cmd.iterate(pairs[0], 'stored.nitro = chain, resi, resn')
                cmd.iterate(pairs[1], 'stored.oxy = chain, resi, resn')
                print ' %s%7s%5d | %s%7s%5d |%7s' % (stored.nitro[0], stored.nitro[2], int(stored.nitro[1]),
                        stored.oxy[0], stored.oxy[2], int(stored.oxy[1]), wrappers)
            sel.append(pairs[0])
            sel.append(pairs[1])
    cmd.delete('_nonpolar')

    if len(sel) > 0:
        cmd.hide('%s' % selection) 
        cmd.show('cartoon','%s' % selection)
        cmd.show_as('dashes', name)
    elif not quiet:
        print ' - no dehydrons were found - '

cmd.extend('dehydron', dehydron)

# vi:expandtab:smarttab
