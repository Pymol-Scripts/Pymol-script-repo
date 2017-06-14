'''
Wrappy: A dehydron calculator plugin for PyMOL
Version: 2.1
Described at PyMOL wiki:
http://www.pymolwiki.org/index.php/dehydron

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: March 2014
License: MIT License
Acknowledgement: The H-bond detection code is based on the list_mc_hbonds.py
script from Robert L. Campbell http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
'''

import sys

if sys.version_info[0] < 3:
    import Tkinter
    from Tkinter import *
else:
    import tkinter as Tkinter
    from tkinter import *

import Pmw
from pymol import cmd
from pymol import stored


def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                             'wrappy',
                             label='Wrappy',
                             command=lambda: mainDialog(self.root))


def mainDialog(root=None):
    """ Creates the GUI"""

    def get_dehydrons():
        angle_range = float(angle_value.get())
        max_distance = float(dist_cutoff_value.get())
        desolv = float(desolv_sphere.get())
        min_wrappers = float(min_value.get())
        max_wrappers = float(max_value.get())
        selection = sel_value.get()
        if len(cmd.get_names()) > 0:
            dehydron(selection, angle_range, max_distance, desolv, min_wrappers, max_wrappers)
        else:
            print('#' * 40)
            print('Please load a protein structure')
            print('#' * 40)
    master = Tkinter.Toplevel(root)
    master.title(' Wrappy ')
    w = Tkinter.Label(master, text='dehydron calculator\nOsvaldo Martin - omarti@unsl.edu.ar',
                      background='#000000',
                      foreground='#cecece',
                      #pady = 20,
                      )
    w.pack(expand=1, fill='both', padx=4, pady=4)

    Pmw.initialise(master)
    nb = Pmw.NoteBook(master, hull_width=420, hull_height=280)
    p1 = nb.add('Main')
    p2 = nb.add('About')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ Main TAB #################################
# hydrogen bond settings
    group = Pmw.Group(p1, tag_text='Hydrogen bond Settings')
    group.pack(fill='x', expand=1, padx=20, pady=1)
    Label(group.interior(), text='angle range').grid(row=2, column=0)
    angle_value = StringVar(master=group.interior())
    angle_value.set(40)
    entry_angle = Entry(group.interior(), textvariable=angle_value, width=10)
    entry_angle.grid(row=2, column=1)
    entry_angle.configure(state='normal')
    entry_angle.update()
    Label(group.interior(), text='max distance').grid(row=3, column=0)
    dist_cutoff_value = StringVar(master=group.interior())
    dist_cutoff_value.set(3.5)
    entry_dist = Entry(group.interior(), textvariable=dist_cutoff_value, width=10)
    entry_dist.grid(row=3, column=1)
    entry_dist.configure(state='normal')
    entry_dist.update()
# dehydron settings
    group = Pmw.Group(p1, tag_text='Dehydron Settings')
    group.pack(fill='x', expand=1, padx=20, pady=5)
    Label(group.interior(), text='desolvatation sphere radius').grid(row=2, column=2)
    desolv_sphere = StringVar(master=group.interior())
    desolv_sphere.set(6.5)
    entry_desolv = Entry(group.interior(), textvariable=desolv_sphere, width=10)
    entry_desolv.grid(row=2, column=3)
    entry_desolv.configure(state='normal')
    entry_desolv.update()
    Label(group.interior(), text='minimum wrappers').grid(row=3, column=2)
    min_value = StringVar(master=group.interior())
    min_value.set(19)
    entry_min_value = Entry(group.interior(), textvariable=min_value, width=10)
    entry_min_value.grid(row=3, column=3)
    entry_min_value.configure(state='normal')
    entry_min_value.update()
    Label(group.interior(), text='maximum wrappers').grid(row=4, column=2)
    max_value = StringVar(master=group.interior())
    max_value.set(35)
    entry_max_value = Entry(group.interior(), textvariable=max_value, width=10)
    entry_max_value.grid(row=4, column=3)
    entry_max_value.configure(state='normal')
    entry_max_value.update()
# selection settings
    group = Pmw.Group(p1, tag_text='Selection')
    group.pack(fill='x', expand=1, padx=20, pady=5)
    Label(group.interior(), text='selection').grid(row=5, column=2)
    sel_value = StringVar(master=group.interior())
    sel_value.set('all')
    entry_sel_value = Entry(group.interior(), textvariable=sel_value, width=10)
    entry_sel_value.grid(row=5, column=3)
    entry_sel_value.configure(state='normal')
    entry_sel_value.update()
# submit
    Button(p1, text="Calculate", command=get_dehydrons).pack(side=BOTTOM)
############################ About TAB #################################
    group = Pmw.Group(p2, tag_text='About wrappy')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    text = u"""For a brief introduction to the dehydron concept, you could
read http://en.wikipedia.org/wiki/dehydron

Citation for this plugin:
Martin O.A.; Wrappy: A dehydron calculator plugin for PyMOL,
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
    text_holder = Text(interior_frame, yscrollcommand=bar.set, foreground="#cecece", background="#000000", font="Times 12")
    bar.config(command=text_holder.yview)
    text_holder.insert(END, text)
    text_holder.pack(side=LEFT, expand="yes", fill="both")
    bar.pack(side=LEFT, expand="yes", fill="y")
    interior_frame.pack(expand="yes", fill="both")


def colorize():
    cmd.hide('(not (name C+CA+N+O))')

    cmd.spectrum('b', 'red_yellow_green', minimum='-1.0', maximum='1.0')
    cmd.select('missing', 'b = -2.0')
    cmd.color('white', 'missing')
    cmd.delete('missing')


def dehydron(selection='all', angle_range=40., max_distance=3.5, desolv=6.5, min_wrappers=19, max_wrappers=35, quiet=0):
    '''
DESCRIPTION

    dehydron calculator

USAGE

    wrappy [ selection [, angle_range [, max_distance [, desolv [, min_wrappers [, max_wrappers ]]]]]]
    '''

    angle, max_distance = float(angle_range), float(max_distance)
    desolv, min_wrappers, max_wrappers = float(desolv), int(min_wrappers), int(max_wrappers)
    quiet = int(quiet)

    DH_name = cmd.get_legal_name('DH_%s' % selection)
    cmd.delete(DH_name)
    HBA_name = cmd.get_legal_name('HBA_%s' % selection)
    cmd.delete(HBA_name)
    HBO_name = cmd.get_legal_name('HBO_%s' % selection)
    cmd.delete(HBO_name)

    selection_hb = '((%s) and polymer)' % (selection)
    hb = cmd.find_pairs("((byres " + selection_hb + ") and n. n)", "((byres " + selection_hb + ") and n. o)", mode=1, cutoff=max_distance, angle=angle_range)

    cmd.select('_nonpolar', '(elem C) and not (solvent or (elem N+O) extend 1)', 0)
    try:
        cmd.select('_selection', '%s' % selection, 0)
    except:
        pass

    low_sel = []
    mean_sel = []
    high_sel = []
    total_wrappers = 0
    for pairs in hb:
        wrappers = cmd.count_atoms('((%s and _nonpolar and _selection) within %f of byca (%s`%d %s`%d))' %
                                   ((pairs[0][0], desolv) + pairs[0] + pairs[1]))
        total_wrappers = total_wrappers + wrappers
        cmd.iterate(pairs[0], 'stored.donor = chain, resi, resn')
        cmd.iterate(pairs[1], 'stored.aceptor = chain, resi, resn')
        if wrappers < min_wrappers:
            cmd.distance(DH_name, pairs[0], pairs[1])
            line = (wrappers, '%12s%4s%6s%6s | %12s%4s%6s%6s |%7s        below' % (pairs[0][0], stored.donor[0], stored.donor[2], stored.donor[1], pairs[1][0], stored.aceptor[0], stored.aceptor[2], stored.aceptor[1], wrappers))
            low_sel.append(line)
        elif wrappers < max_wrappers:
            cmd.distance(HBA_name, pairs[0], pairs[1])
            line = (wrappers, '%12s%4s%6s%6s | %12s%4s%6s%6s |%7s      average' % (pairs[0][0], stored.donor[0], stored.donor[2], stored.donor[1], pairs[1][0], stored.aceptor[0], stored.aceptor[2], stored.aceptor[1], wrappers))
            mean_sel.append(line)
        else:
            cmd.distance(HBO_name, pairs[0], pairs[1])
            line = (wrappers, '%12s%4s%6s%6s | %12s%4s%6s%6s |%7s         over' % (pairs[0][0], stored.donor[0], stored.donor[2], stored.donor[1], pairs[1][0], stored.aceptor[0], stored.aceptor[2], stored.aceptor[1], wrappers))
            high_sel.append(line)

    cmd.delete('_nonpolar')
    cmd.delete('_selection')
    # compute the z_scores. Useful for protein structure validation.
    stored.ResiduesNames = []
    cmd.iterate('(name ca)', 'stored.ResiduesNames.append((resn))')
    total_residues = float(len(stored.ResiduesNames))
    z_score_wrappers = ((total_wrappers / total_residues) - 17) / 2
    z_score_hb = ((len(hb) / total_residues) - 0.62) / 0.06

    if len(low_sel) > 0:
        cmd.show_as('dashes', DH_name)
        cmd.color('red', DH_name)
        low_sel.sort()
    if len(mean_sel) > 0:
        cmd.show_as('dashes', HBA_name)
        cmd.color('yellow', HBA_name)
        mean_sel.sort()
    if len(high_sel) > 0:
        cmd.show_as('dashes', HBO_name)
        cmd.color('green', HBO_name)
        high_sel.sort()

    if not quiet:
        hb.sort()
        print("--------------------------------------------------------------------")
        print("------------------------------Results ------------------------------")
        print("--------------------------------------------------------------------")
        print("           Donor             |            Aceptor           |")
        print("    Object   Chain Residue   |     Object   Chain Residue   | # wrappers  wrapping")
        for line in low_sel:
            print(line[1])
        for line in mean_sel:
            print(line[1])
        for line in high_sel:
            print(line[1])
        print('\nProtein global statistics:')
        print('\nz-score wrappers = %6.2f\nz-score hydrogen bonds = %6.2f\n' % (z_score_wrappers, z_score_hb))
    elif not quiet and len(hb) != 0:
        print('\n - no dehydrons were found - ')
        print('\nz-score hydrogen bonds = %6.2f\n' % (z_score_hb))
    else:
        print('\n - no hydrogen bonds were found - ')

cmd.extend('wrappy', dehydron)

# vi:expandtab:smarttab
