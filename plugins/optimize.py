'''
Optimize
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/optimize

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: October 2013
License: GNU General Public License
Version 0.1
'''

import Tkinter
from Tkinter import *
import Pmw
from pymol import cmd
import openbabel as ob


def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'Optimize',
                            label = 'Optimize',
                            command = lambda : mainDialog())

def mainDialog():
    """ Creates the GUI """
    global entry_vdw, entry_elec
    def set_minimize():
        forcefield = ff_value.get()
        method = method_value.get()
        nsteps = int(entry_nsteps.get())
        conv = float(entry_conv.get())
        cutoff = bool(cutoff_value.get())
        cut_vdw = float(entry_vdw.get())
        cut_elec = float(entry_elec.get())
        selection = sel_value.get()
        minimize(selection, forcefield, method, nsteps, conv, cutoff, cut_vdw, cut_elec)

    master = Tk()
    master.title(' Optimize ')
    w = Tkinter.Label(master, text="\nOptimize: Let's find that minima!\n",
                                background = 'black',
                                foreground = 'white')
    w.pack(expand=1, fill = 'both', padx=4, pady=4)
############################ NoteBook #########################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=430, hull_height=320)
    p1 = nb.add(' Local optimization ')
    p2 = nb.add(' Global Optimization ')
    p3 = nb.add('    About   ')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ Minimization TAB #################################
    group = Pmw.Group(p1,tag_text='Minimization options')
    group.pack(fill='both', expand=1, padx=5, pady=5)
# Force Field options
    ff_value = StringVar(master=group.interior())
    ff_value.set('GAFF')
    ff_menu = Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = 'Force Field',
                menubutton_textvariable = ff_value,
                items = ['GAFF', 'MMFF94s', 'MMFF94', 'UFF', 'Ghemical'],
                menubutton_width = 15,
        ).grid(row=0, columnspan=2)
# Method
    method_value = StringVar(master=group.interior())
    method_value.set('conjugate gradients')
    Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = '   Method  ',
                menubutton_textvariable = method_value,
                items = ['conjugate gradients', 'steepest descent'],
                menubutton_width = 15,
        ).grid(row=1, columnspan=2)
#
    Label(group.interior(), text='steps').grid(row=2, column=0)
    nsteps = StringVar(master=group.interior())
    nsteps.set(300)
    entry_nsteps = Entry(group.interior(),textvariable=nsteps, width=15)
    entry_nsteps.grid(row=2, column=1)
    entry_nsteps.update()
    Label(group.interior(), text='convergence').grid(row=3, column=0)
    conv = StringVar(master=group.interior())
    conv.set(0.01)
    entry_conv = Entry(group.interior(),textvariable=conv, width=15)
    entry_conv.grid(row=3, column=1)
    entry_conv.update()
    Label(group.interior(), text='selection').grid(row=4, column=0)
    sel_value = StringVar(master=group.interior())
    names = cmd.get_names('all')
    if len(names) > 0:
        sel_value.set(names[0])
    else:
        sel_value.set('')
    entry_sel_value = Entry(group.interior(),textvariable=sel_value, width=15)
    entry_sel_value.grid(row=4, column=1)
    entry_sel_value.configure(state='normal')
    entry_sel_value.update()

###########################################################################
###########################################################################
###########################################################################
    cutoff_value = BooleanVar(master=group.interior())
    cutoff_value.set(False)
    Radiobutton(group.interior(), text='No cutoff ', variable=cutoff_value, value=False, command=disable_entry).grid(row=5, columnspan=3)
    Radiobutton(group.interior(), text='Use cutoff', variable=cutoff_value, value=True,
command=enable_entry).grid(row=6, columnspan=3)
# http proxy and Port
    Label(group.interior(), text='Van der Waals').grid(row=7, column=0)
    vdw_value = StringVar(master=group.interior())
    vdw_value.set(6.0)
    entry_elec = Entry(group.interior(),textvariable=vdw_value, width=15)
    entry_elec.grid(row=7, column=1)
    entry_elec.configure(state='disabled')
    entry_elec.update()
    Label(group.interior(), text='Electrostatic').grid(row=8, column=0)
    elec_value = StringVar(master=group.interior())
    elec_value.set(8.0)
    entry_vdw = Entry(group.interior(),textvariable=elec_value, width=15)
    entry_vdw.grid(row=8, column=1)
    entry_vdw.configure(state='disabled')
    entry_vdw.update()

# Run
    Button(p1, text="Minimize", command=set_minimize).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text =u"""
Monte Carlo Routine (coming soon)
""",justify=LEFT).pack()
#    Button(p2, text="Reset View", command=colorize).pack(side=BOTTOM)
############################ About TAB ########################################
    Label(p3, text = """
Optimize provides a PyMOL graphical interface to some 
of the many options available in openbabel (openbabel.org).


If you find Optimize useful great! 
If you don't and have some suggestions or coments 
to doplease write to me (aloctavodia@gmail.com).
""",justify=CENTER).pack()
    master.mainloop()


def enable_entry():
    """enables the fields for proxy and port"""
    entry_vdw.configure(state='normal')
    entry_vdw.update()
    entry_elec.configure(state='normal')
    entry_elec.update()


def disable_entry():
    """disables all the fields related to the proxy tab"""
    entry_vdw.configure(state='disabled')
    entry_vdw.update()
    entry_elec.configure(state='disabled')
    entry_elec.update()


def minimize(selection='all', forcefield='GAFF', method='conjugate gradients', nsteps= 300, conv=0.01, cutoff=False, cut_vdw=6.0, cut_elec=8.0):
    cmd.h_add(selection)
    pdb_string = cmd.get_pdbstr(selection)
    name = cmd.get_legal_name(selection)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('pdb', 'pdb')
    mol = ob.OBMol()
    obconversion.ReadString(mol, pdb_string)
    ff = ob.OBForceField.FindForceField(forcefield) ## GAFF, MMFF94s, MMFF94, UFF, Ghemical
    ff.Setup(mol)
    if cutoff == True:
        ff.EnableCutOff(True)
        ff.SetVDWCutOff(cut_vdw)
        ff.SetElectrostaticCutOff(cut_elec)
    if method == 'conjugate gradients':
        ff.ConjugateGradients(nsteps, conv)
    else:
        ff.SteepestDescent(nsteps, conv)
    ff.GetCoordinates(mol)
    nrg = ff.Energy()
    pdb_string = obconversion.WriteString(mol)
    cmd.delete(name)
    if name == 'all':
        name = 'all_'
    cmd.read_pdbstr(pdb_string, name)
#    print '#########################################'
#    print 'The Energy of %s is %8.2f kcal/mol       '  % (selection, nrg)
#    print '#########################################'

cmd.extend('minimize', minimize)
