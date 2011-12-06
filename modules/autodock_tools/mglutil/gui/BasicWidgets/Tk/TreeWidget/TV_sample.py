from TreeWidget.tree import TreeView

if __name__ == '__main__':
    tv = TreeView()
    # addNode(nodename, parentname = None)
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
   

    tv.addNode('protein_2')
    tv.addNode('protein_3')


    tv.addNode('residue_21',    parent='protein_2')
    tv.addNode('residue_25',    parent='protein_2')
    tv.addNode('basdfe',        parent='protein_2|residue_21')
    tv.addNode('AminoAcid',     parent='protein_2|residue_21')
   
 
    tv.addNode('etc',       parent='protein_1|residue_11')
    tv.addNode('color',     parent='protein_1|residue_11|etc')
    tv.addNode('density',   parent='protein_1|residue_11|etc')
    tv.addNode('residue_12',parent='protein_1')
    
  
    for a in range(1):
        name = 'AA' + str(a)
        tv.addNode(name, parent='protein_1|residue_11|AminoAcid')


    tv.addNode('2', parent='protein_2|residue_21')
    tv.addNode('3', parent='protein_2|residue_21')
    tv.addNode('4', parent='protein_2|residue_21')

    tv.addNode('L', parent='protein_2|residue_21|AminoAcid')
    tv.addNode('S', parent='protein_2|residue_21|AminoAcid')
 
    
    for a in range(10):
        name = 'A' + str(a)
        tv.addNode(name,  parent='protein_2|residue_21|AminoAcid')


    tv.addNode('protein_4')
    tv.addNode('residue_22', parent='protein_2')

##  to delete a node:
#   tv.deleteNode(nodename, parentname)
# e.g. >>> tv.deleteNode('residue_21', 'protein_2')
# e.g. >>> tv.deleteNode('protein_2')
