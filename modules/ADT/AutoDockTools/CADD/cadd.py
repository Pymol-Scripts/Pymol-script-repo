import os

class CADD_Vision:
    """
    This class supports creating Vision network using Opal web services to
    perform computational vistual screening experiments
    """

    def __init__(self):

        # create an instance of Vision without a GUI
        from Vision.VPE import NoGuiExec
        self.ed = NoGuiExec()

    def run(self):
        from NetworkEditor.net import Communicator
        self.net.communicator = Communicator(self.net)
        print 'Communicator listening on port:', self.net.communicator.port

        import socket
        #f = open(argv[0]+'.sock', 'w')
        #f.write("%s %i"%(socket.gethostbyname(socket.gethostname()),
        #                 self.net.communicator.port))
        #f.close()
        print "%s %i"%(socket.gethostbyname(socket.gethostname()),
                       self.net.communicator.port)
                       
        self.net.run()


class CADD_VistualScreening(CADD_Vision):
    """
    virtual screenign based on VS_net.py

    This network allows a user to specify:
      - a directory of structures
      - a set of ligands taken from a public ligand database
      - filtering options for the ligands
      - whether to preserve the charges or not

    the parameters are set using hte correspondign methods
        .setStructureDirectory(path)
        .setLigandDB(string)
        .setFilter(string)
        .setPreseveCharges(0/1)
    """

    def __init__(self):

        CADD_Vision.__init__(self)

        # name used by Vision for this network
        name = 'vs_opal'
        
        # find the path to CADD networks in the AutoDockTools package
        from mglutil.util.packageFilePath import findFilePath
        netpath = findFilePath('CADD', 'AutoDockTools')

        # build full path to VS_net.py network
        filename = os.path.join(netpath, 'VisionNetworks', 'VS_net.py')

        # load the network
        self.net = self.ed.loadNetworkNOGUI(name, filename)


    def setStructureDirectory(self, path):
        assert os.path.exists(path)
        assert os.path.isdir(path)

        # get handle to input ports
        nodes = self.net.getNodeByName('GetStructuresFromDir')
        assert len(nodes)==1
        node = nodes[0]

        # a get a handle to the port
        port = node.getInputPortByName('directory')

        # set the widget to the specified value without triggering execution
        port.widget.set(path, run=0)

        
    def setLigandDB(self, name):
        assert name in [
            'sample', 'sample2', 'NCI_DS1', 'NCI_DS2', 'NCIDS_SC', 'oldNCI',
            'chembridge_building_blocks', 'drugbank_neutraceutics',
            'drugbank_smallmol', 'fda_approved', 'human_metabolome',
            'asinex', 'otava', 'zinc_natural_products.prop']

        # get handle to input ports
        nodes = self.net.getNodeByName('PublicServerLigandDB')
        assert len(nodes)==1
        node = nodes[0]

        # a get a handle to the port
        port = node.getInputPortByName('server_lib')

        # set the widget to the specified value without triggering execution
        port.widget.set(name, run=0)


    def setFilter(self, filterMode):
        assert filterMode in ['None', 'default', 'Lipinski-like',
                              'Drug-like', 'Drug-like frag', 'Cistom']
        # get handle to input ports
        nodes = self.net.getNodeByName('FilterLigandsNode')
        assert len(nodes)==1
        node = nodes[0]

        # a get a handle to the port
        port = node.getInputPortByName('filterMode')

        # set the widget to the specified value without triggering execution
        port.widget.set(filterMode, run=0)


    def setPreseveCharges(self, onOff):
        assert onOff in (0,1)
        
        # get handle to input ports
        nodes = self.net.getNodeByName('PreserveCharges?')
        assert len(nodes)==1
        node = nodes[0]

        # a get a handle to the port
        port = node.getInputPortByName('button')

        # set the widget to the specified value without triggering execution
        port.widget.set(onOff, run=0)


