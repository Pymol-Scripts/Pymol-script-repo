#########################################################################
#
# Date: Aug 2003  Authors: Daniel Stoffler, Michel Sanner
#
#       stoffler@scripps.edu
#       sanner@scripps.edu
#
# Copyright: Daniel Stoffler, Michel Sanner, and TSRI
#
#########################################################################

import sys
ed = None

def setUp():
    global ed
    from Vision.VPE import VisualProgramingEnvironment
    ed = VisualProgramingEnvironment(name='Vision', withShell=0,)
    ed.master.update_idletasks()
    ed.configure(withThreads=0)

def tearDown():
    ed.exit_cb()
    import gc
    gc.collect()

##########################
## Helper methods
##########################

def pause(sleepTime=0.4):
    from time import sleep
    ed.master.update()
    sleep(sleepTime)

##########################
## Tests
##########################

def test_01_loadMolKitLib():
    from MolKit.VisionInterface.MolKitNodes import molkitlib
    ed.addLibraryInstance(molkitlib, 'MolKit.VisionInterface.MolKitNodes', 'molkitlib')
    ed.root.update_idletasks()
    pause()


def test_02_allMolKitNodes():
    from MolKit.VisionInterface.MolKitNodes import molkitlib
    ed.addLibraryInstance(
        molkitlib, 'MolKit.VisionInterface.MolKitNodes', 'molkitlib')
    ed.root.update_idletasks()
    pause()
    # test the molkit nodes
    lib = 'MolKit'
    libs = ed.libraries
    posx = 150
    posy = 150

    #ed.ModulePages.selectpage(lib)
    ed.root.update_idletasks()
    for cat in libs[lib].libraryDescr.keys():
        for node in libs[lib].libraryDescr[cat]['nodes']:
            klass = node.nodeClass
            kw = node.kw
            args = node.args
            netNode = apply( klass, args, kw )
            print 'testing: '+node.name # begin node test
            #add node to canvas
            ed.currentNetwork.addNode(netNode,posx,posy)
            # show widget in node if available:
            widgetsInNode = netNode.getWidgetsForMaster('Node')
            if len( widgetsInNode.items() ):
                if not netNode.isExpanded():
                    netNode.toggleNodeExpand_cb()
                    ed.root.update_idletasks()
                # and then hide it
                netNode.toggleNodeExpand_cb()
                ed.root.update_idletasks()

            # show widgets in param panel if available:
            widgetsInPanel = netNode.getWidgetsForMaster('ParamPanel')
            if len(widgetsInPanel.items()):
                netNode.paramPanel.show()
                ed.root.update_idletasks()
                #and then hide it
                netNode.paramPanel.hide()
                ed.root.update_idletasks()

            # and now delete the node
            ed.currentNetwork.deleteNodes([netNode])
            ed.root.update_idletasks()
            print 'passed: '+node.name # end node test


def test_03_cpkNet():
    ed.loadNetwork('cpk_net.py')
    pause()

    ed.runCurrentNet_cb()
    pause()
    
