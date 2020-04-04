#!/usr/bin/python

'''PyMOL plugin that provides show_contacts command and GUI 
for highlighting good and bad polar contacts. Factored out of 
clustermols by Matthew Baumgartner.
The advantage of this package is it requires many fewer dependencies.
'''
from __future__ import print_function

import sys,os
from pymol import cmd

DEBUG=1

def show_contacts(selection,selection2,result="contacts",cutoff=3.6, bigcutoff = 4.0, SC_DEBUG = DEBUG):
    """
    USAGE
    
    show_contacts selection, selection2, [result=contacts],[cutoff=3.6],[bigcutoff=4.0]
    
    Show various polar contacts, the good, the bad, and the ugly.
    
    Edit MPB 6-26-14: The distances are heavy atom distances, so I upped the default cutoff to 4.0
    
    Returns:
    True/False -  if False, something went wrong
    """
    if SC_DEBUG > 4:
        print('Starting show_contacts')
        print('selection = "' + selection + '"')
        print('selection2 = "' + selection2 + '"')
            
    result = cmd.get_legal_name(result)

    #if the group of contacts already exist, delete them
    cmd.delete(result)

    # ensure only N and O atoms are in the selection
    all_don_acc1 = selection + " and (donor or acceptor)"
    all_don_acc2 = selection2 + " and  (donor or acceptor)"
    
    if SC_DEBUG > 4:
        print('all_don_acc1 = "' + all_don_acc1 + '"')
        print('all_don_acc2 = "' + all_don_acc2 + '"')
    
    #if theses selections turn out not to have any atoms in them, pymol throws cryptic errors when calling the dist function like:
    #'Selector-Error: Invalid selection name'
    #So for each one, manually perform the selection and then pass the reference to the distance command and at the end, clean up the selections
    #the return values are the count of the number of atoms
    all1_sele_count = cmd.select('all_don_acc1_sele', all_don_acc1)
    all2_sele_count = cmd.select('all_don_acc2_sele', all_don_acc2)
    
    #print out some warnings
    if DEBUG > 3:
        if not all1_sele_count:
            print('Warning: all_don_acc1 selection empty!')
        if not all2_sele_count:
            print('Warning: all_don_acc2 selection empty!')
    
    ########################################
    allres = result + "_all"
    if all1_sele_count and all2_sele_count:
        cmd.distance(allres, 'all_don_acc1_sele', 'all_don_acc2_sele', bigcutoff, mode = 0)
        cmd.set("dash_radius", "0.05", allres)
        cmd.set("dash_color", "purple", allres)
        cmd.hide("labels", allres)
    
    ########################################
    #compute good polar interactions according to pymol
    polres = result + "_polar"
    if all1_sele_count and all2_sele_count:
        cmd.distance(polres, 'all_don_acc1_sele', 'all_don_acc2_sele', cutoff, mode = 2) #hopefully this checks angles? Yes
        cmd.set("dash_radius","0.126",polres)
    
    ########################################
    #When running distance in mode=2, the cutoff parameter is ignored if set higher then the default of 3.6
    #so set it to the passed in cutoff and change it back when you are done.
    old_h_bond_cutoff_center = cmd.get('h_bond_cutoff_center') # ideal geometry
    old_h_bond_cutoff_edge = cmd.get('h_bond_cutoff_edge') # minimally acceptable geometry
    cmd.set('h_bond_cutoff_center', bigcutoff)
    cmd.set('h_bond_cutoff_edge', bigcutoff)
        
    #compute possibly suboptimal polar interactions using the user specified distance
    pol_ok_res = result + "_polar_ok"
    if all1_sele_count and all2_sele_count:
        cmd.distance(pol_ok_res, 'all_don_acc1_sele', 'all_don_acc2_sele', bigcutoff, mode = 2) 
        cmd.set("dash_radius", "0.06", pol_ok_res)

    #now reset the h_bond cutoffs
    cmd.set('h_bond_cutoff_center', old_h_bond_cutoff_center)
    cmd.set('h_bond_cutoff_edge', old_h_bond_cutoff_edge) 
    
    
    ########################################
    
    onlyacceptors1 = selection + " and (acceptor and !donor)"
    onlyacceptors2 = selection2 + " and (acceptor and !donor)"
    onlydonors1 = selection + " and (!acceptor and donor)"
    onlydonors2 = selection2 + " and (!acceptor and donor)"  
    
    #perform the selections
    onlyacceptors1_sele_count = cmd.select('onlyacceptors1_sele', onlyacceptors1)
    onlyacceptors2_sele_count = cmd.select('onlyacceptors2_sele', onlyacceptors2)
    onlydonors1_sele_count = cmd.select('onlydonors1_sele', onlydonors1)
    onlydonors2_sele_count = cmd.select('onlydonors2_sele', onlydonors2)    
    
    #print out some warnings
    if SC_DEBUG > 2:
        if not onlyacceptors1_sele_count:
            print('Warning: onlyacceptors1 selection empty!')
        if not onlyacceptors2_sele_count:
            print('Warning: onlyacceptors2 selection empty!')
        if not onlydonors1_sele_count:
            print('Warning: onlydonors1 selection empty!')
        if not onlydonors2_sele_count:
            print('Warning: onlydonors2 selection empty!')    
            
    
    accres = result+"_aa"
    if onlyacceptors1_sele_count and onlyacceptors2_sele_count:
        aa_dist_out = cmd.distance(accres, 'onlyacceptors1_sele', 'onlyacceptors2_sele', cutoff, 0)

        if aa_dist_out < 0:
            print('\n\nCaught a pymol selection error in acceptor-acceptor selection of show_contacts')
            print('accres:', accres)
            print('onlyacceptors1', onlyacceptors1)
            print('onlyacceptors2', onlyacceptors2)
            return False
    
        cmd.set("dash_color","red",accres)
        cmd.set("dash_radius","0.125",accres)
    
    ########################################
    
    donres = result+"_dd"
    if onlydonors1_sele_count and onlydonors2_sele_count:
        dd_dist_out = cmd.distance(donres, 'onlydonors1_sele', 'onlydonors2_sele', cutoff, 0)
        
        #try to catch the error state 
        if dd_dist_out < 0:
            print('\n\nCaught a pymol selection error in dd selection of show_contacts')
            print('donres:', donres)
            print('onlydonors1', onlydonors1)
            print('onlydonors2', onlydonors2)
            print("cmd.distance('" + donres + "', '" + onlydonors1 + "', '" + onlydonors2 + "', " + str(cutoff) + ", 0)")  
            return False
        
        cmd.set("dash_color","red",donres)  
        cmd.set("dash_radius","0.125",donres)
    
    ##########################################################
    ##### find the buried unpaired atoms of the receptor #####
    ##########################################################
    
    #initialize the variable for when CALC_SASA is False
    unpaired_atoms = ''
    
        
    ## Group
    cmd.group(result,"%s %s %s %s %s %s" % (polres, allres, accres, donres, pol_ok_res, unpaired_atoms))
    
    ## Clean up the selection objects
    #if the show_contacts debug level is high enough, don't delete them.
    if SC_DEBUG < 5:
        cmd.delete('all_don_acc1_sele')
        cmd.delete('all_don_acc2_sele')
        cmd.delete('onlyacceptors1_sele')
        cmd.delete('onlyacceptors2_sele')
        cmd.delete('onlydonors1_sele')
        cmd.delete('onlydonors2_sele')
    
    
    return True
cmd.extend('contacts', show_contacts) #contacts to avoid clashing with cluster_mols version



    
    
#################################################################################
########################### Start of pymol plugin code ##########################
#################################################################################


about_text = '''show_contacts was factored out of the much more full-featured cluster_mols
by Dr. Matt Baumgartner (https://pymolwiki.org/index.php/Cluster_mols).  It provides
an easy way to highlight polar contacts (and clashes) between two selections without
requiring the installation of additional dependencies.
'''

class Show_Contacts:
    ''' Tk version of the Plugin GUI '''
    def __init__(self, app):
        parent = app.root
        self.parent = parent
        
        self.app = app
        
        import Pmw

        ############################################################################################
        ### Open a window with options to select to loaded objects ###
        ############################################################################################

        self.select_dialog = Pmw.Dialog(parent, 
                         buttons = ('Ok','Cancel'), 
                         title = 'Show Contacts Plugin',
                         command = self.button_pressed )
    
        self.select_dialog.withdraw()
    

        #allow the user to select from objects already loaded in pymol
        self.select_object_combo_box = Pmw.ComboBox(self.select_dialog.interior(),
                                                               scrolledlist_items=[],
                                                               labelpos='w',
                                                               label_text='Select loaded object:',
                                                               listbox_height = 2,
                                                               dropdown=True)
        self.select_object_combo_box2 = Pmw.ComboBox(self.select_dialog.interior(),
                                                               scrolledlist_items=[],
                                                               labelpos='w',
                                                               label_text='Select loaded object:',
                                                               listbox_height = 2,
                                                               dropdown=True)                                                               
        self.select_object_combo_box.grid(column=1, row=0)
        self.select_object_combo_box2.grid(column=2, row=0)
        self.populate_ligand_select_list()
        self.select_dialog.show()
        

      
    def button_pressed(self, result):
        if hasattr(result,'keycode'):
            if result.keycode == 36:
                print('keycode:', result.keycode)
        elif result == 'Ok' or result == 'Exit' or result == None:
            s1 = self.select_object_combo_box.get()
            s2 = self.select_object_combo_box2.get()
            show_contacts(s1,s2,'%s_%s'%(s1,s2))
            self.select_dialog.withdraw()            
        elif result == 'Cancel' or result == None:
            self.select_dialog.withdraw()

            

    
    def populate_ligand_select_list(self):
        ''' Go thourgh the loaded objects in PyMOL and add them to the selected list. '''
        #get the loaded objects
        loaded_objects = _get_select_list()
         
        self.select_object_combo_box.clear()
        self.select_object_combo_box2.clear()
        
        for ob in loaded_objects:
            self.select_object_combo_box.insert('end', ob)
            self.select_object_combo_box2.insert('end', ob)
        

def _get_select_list():
    '''
    Get either a list of object names, or a list of chain selections
    '''
    loaded_objects = [name for name in cmd.get_names('all', 1) if '_cluster_' not in name]

    # if single object, try chain selections
    if len(loaded_objects) == 1:
        chains = cmd.get_chains(loaded_objects[0])
        if len(chains) > 1:
            loaded_objects = ['{} & chain {}'.format(loaded_objects[0], chain) for chain in chains]

    return loaded_objects


class Show_Contacts_Qt_Dialog(object):
    ''' Qt version of the Plugin GUI '''
    def __init__(self):
        from pymol.Qt import QtWidgets
        dialog = QtWidgets.QDialog()
        self.setupUi(dialog)
        self.populate_ligand_select_list()
        dialog.accepted.connect(self.accept)
        dialog.exec_()

    def accept(self):
        s1 = self.select_object_combo_box.currentText()
        s2 = self.select_object_combo_box2.currentText()
        show_contacts(s1, s2, '%s_%s' % (s1, s2))

    def populate_ligand_select_list(self):
        loaded_objects = _get_select_list()

        self.select_object_combo_box.clear()
        self.select_object_combo_box2.clear()

        self.select_object_combo_box.addItems(loaded_objects)
        self.select_object_combo_box2.addItems(loaded_objects)

        if len(loaded_objects) > 1:
            self.select_object_combo_box2.setCurrentIndex(1)

    def setupUi(self, Dialog):
        # Based on auto-generated code from ui file
        from pymol.Qt import QtCore, QtWidgets
        Dialog.resize(400, 50)
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        label = QtWidgets.QLabel("Select loaded object:", Dialog)
        self.gridLayout.addWidget(label, 0, 0, 1, 1)
        self.select_object_combo_box = QtWidgets.QComboBox(Dialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.select_object_combo_box.setSizePolicy(sizePolicy)
        self.select_object_combo_box.setEditable(True)
        self.gridLayout.addWidget(self.select_object_combo_box, 0, 1, 1, 1)
        label = QtWidgets.QLabel("Select loaded object:", Dialog)
        self.gridLayout.addWidget(label, 1, 0, 1, 1)
        self.select_object_combo_box2 = QtWidgets.QComboBox(Dialog)
        self.select_object_combo_box2.setSizePolicy(sizePolicy)
        self.select_object_combo_box2.setEditable(True)
        self.gridLayout.addWidget(self.select_object_combo_box2, 1, 1, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.gridLayout.addWidget(self.buttonBox, 2, 0, 1, 2)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)

    
def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Show Contacts', Show_Contacts_Qt_Dialog)
        return
    except Exception as e:
        print(e)
    self.menuBar.addmenuitem('Plugin', 'command', 'Show Contacts', label = 'Show Contacts', command = lambda s=self : Show_Contacts(s))  
        

