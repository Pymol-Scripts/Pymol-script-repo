##############################################################################
# dehydron:  A dehydron calculator plugin for PyMOL
# Author: Osvaldo Martin
# e-mail: aloctavodia@gmail.com
# License: GNU General Public License
#
# Acknowledgement: The H-bond detection code is based on the list_mc_hbonds.py 
# script from Robert L. Campbell http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
#
##############################################################################

import Tkinter
from Tkinter import *
import Pmw
from pymol import cmd
from pymol import stored # Import PyMOL's stored module.  This will allow us with a 
# way to pull out the PyMOL data and modify it in our script.
# See below.


def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'dehydron',
                            label = 'dehydron',
                            command = lambda : mainDialog())


def mainDialog():
    """ Creates the GUI"""
    global angle_value
    global dist_cutoff_value
    global desolv_sphere
    global min_value
    master = Tkinter.Tk()
    master.title(' dehydron ')
    w = Tkinter.Label(master, text = 'Dehydron calculator\nOsvaldo Martin - omarti@unsl.edu.ar',
                                background = '#000000',
                                foreground = '#cecece',
                                #pady = 20,
                                )
    w.pack(expand=1, fill = 'both', padx = 4, pady = 4)

    Pmw.initialise(master)
    nb = Pmw.NoteBook(master, hull_width = 420, hull_height=240)
    p1 = nb.add('Main')
    p2 = nb.add('About')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################  Main TAB  #################################
# hydrogen bond settings
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
#dehydron settings
    group = Pmw.Group(p1,tag_text='Dehydron Settings')
    group.pack(fill='x', expand=1, padx=20, pady=5)
    Label(group.interior(), text='desolvatation sphere radius').grid(row=2, column=2)
    desolv_sphere = StringVar(master=group.interior())
    desolv_sphere.set(6.5)
    entry_desolv=Entry(group.interior(),textvariable=desolv_sphere, width=10)
    entry_desolv.grid(row=2, column=3)
    entry_desolv.configure(state='normal')
    entry_desolv.update()
    Label(group.interior(), text='cut-off dehydrons').grid(row=3, column=2)
    min_value = StringVar(master=group.interior())
    min_value.set(19)
    entry_min_value=Entry(group.interior(),textvariable=min_value, width=10)
    entry_min_value.grid(row=3, column=3)
    entry_min_value.configure(state='normal')
    entry_min_value.update()
# submit
    Button(p1, text="Calculate", command=get_dehydrons).pack(side=BOTTOM)
############################ About TAB #################################
    group = Pmw.Group(p2, tag_text='About dehydron plug-in')
    group.pack(fill = 'both', expand=1, padx = 5, pady = 5)
    text =u"""For a brief introduction to the dehydron concept, you could
read http://en.wikipedia.org/wiki/Dehydron

Citation for this plugin:
    Martin O.A.; Dehydron calculator plugin for PyMOL, 
2012. IMASL-CONICET.

Citation for PyMOL may be found here:
    http://pymol.sourceforge.net/faq.html#CITE

Citation for Dehydrons (I think these could be used):
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
    # Add this as text in a scrollable pane.
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


def get_dehydrons():
    cmd.delete('dehydrons')
    cmd.delete('DH_pairs')
    cmd.hide()
    angle = float(angle_value.get())
    cutoff = float(dist_cutoff_value.get())
    desolv = float(desolv_sphere.get())
    min_wrappers = float(min_value.get())
    selection = 'name n or name o and not resn hoh'
    hb = cmd.find_pairs("((byres "+selection+") and n. n)","((byres "+selection+") and n. o)",mode=1,cutoff=cutoff,angle=angle)
# sort the list for easier reading
    hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))
    print "------------------------------------------------\n----------------dehydron Results-----------------\n------------------------------------------------\n      Donor     |      Aceptor      | \nChain   Residue | Chain   Residue   | # dehydrons"
    sel = []
    wra = 0
    for pairs in hb:
        cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]), 'stored.nitro = chain, resi, resn') # extracts the nitrogen
        cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]), 'stored.oxy = chain, resi, resn') # extracts the oxygen
        wrappers = cmd.select('wrap', '(((chain %s and name ca and resi %s) around %f) or ((chain %s and name ca and resi %s) around %f)) and (not ((neighbor name n*+o*) or (name n*+o*+h*)))' % (stored.nitro[0], stored.nitro[1], desolv, stored.oxy[0], stored.oxy[1], desolv)) #Identifies the putative dehydrons
        if wrappers < min_wrappers:
            wra = 1
            cmd.distance('Dehydrons',"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
            print '  %s%7s%5d |   %s%7s%5d   |%7s' % (stored.nitro[0], stored.nitro[2], int(stored.nitro[1]), stored.oxy[0], stored.oxy[2], int(stored.oxy[1]), wrappers)
            if stored.nitro[1] not in sel:
                sel.append(stored.nitro[1])
            if stored.oxy[1] not in sel:
                sel.append(stored.oxy[1])
    if wra == 1:
        cmd.select('DH_pairs', 'resi %s' % ('+').join(sel))
        cmd.show('lines', 'DH_pairs')
    cmd.disable('DH_pairs')
    cmd.hide('labels')
#    cmd.delete('wrap')
    cmd.show('cartoon')
    cmd.show('dashes')


#cmd.extend('dehydron', get_dehydrons)
