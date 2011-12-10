from NetworkEditor.items import NetworkNode
class GetStructuresFromDir(NetworkNode):
    """
    Get list of structure file names from directory and create
    objects based on these files
    """
    mRequiredTypes = {}
    mRequiredSynonyms = [
    ]

    def __init__(self, constrkw = {},  name='GetStructuresFromDir', **kw):
        kw['constrkw'] = constrkw
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr
        ip.append({'name': 'directory', 'datatype': 'string'})

        self.widgetDescr['directory'] = {
            'class':'NEEntryWithDirectoryBrowser', 'master':'node', 'width':20,
            'initialValue':'', 'labelCfg':{'text':'Structure Directory: '}
            }

        op = self.outputPortsDescr
        op.append({'name': 'structure_list_obj', 'datatype': 'list'})

#        op = self.outputPortsDescr
#        op.append({'name': 'structure_list', 'datatype': 'list'})


        code = """def doit(self, directory):
        import glob, os
        from AutoDockTools.VisionInterface.Adt.receptor import receptor

        directory = os.path.abspath(directory)

        if directory == None:
            return 'stop'

        if directory is not None:
            cwd = os.getcwd()
            os.chdir(directory)
        try:
            filenames_pdbqt = glob.glob('*.pdbqt')
            filenames_pqr = glob.glob('*.pqr')
            filenames_pdb = glob.glob('*.pdb')

            structure_list = filenames_pdbqt

            for f in filenames_pqr:
                sid = f.rstrip('.pqr')
                s = sid + '.pdbqt'
                found = False
 
                for i in structure_list:
                    if i == s:
                        found = True        
                        break

                if found == False:
                    structure_list.append(f)

            for f in filenames_pdb:
                sid = f.rstrip('.pdb')
                found = False
 
                for i in structure_list:
                    if i == sid + '.pdb' or i == sid + '.pdbqt':
                        found = True        
                        break

                if found == False:
                    structure_list.append(f)

            structure_list = [os.path.join(directory, x) for x in filenames_pdbqt]
            structure_list_obj = [receptor(x) for x in structure_list]
        finally:
            if directory is not None:
                os.chdir(cwd)

        print "--------------------------------------------------------------"
        print "The following structures are found in the structure directory:"

        count = 0

        for i in structure_list:
            count = count + 1
            print str(count) + '. ' + i
           
        print "--------------------------------------------------------------"
        
        pass
#        self.outputData(structure_list_obj=structure_list_obj, structure_list=structure_list)
        self.outputData(structure_list_obj=structure_list_obj)
"""
        self.configure(function=code)


    def beforeAddingToNetwork(self, net):
        try:
            ed = net.getEditor()
        except:
            import traceback; traceback.print_exc()
            print 'Warning! Could not import widgets'

