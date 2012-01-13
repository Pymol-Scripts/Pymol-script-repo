'''
See more at: http://www.pymolwiki.org/index.php/annocryst

######################################################
#
#  AnnoCryst for PyMOL
#
#  Anna Gerber
#  email: agerber@itee.uq.edu.au
#
#  Copyright 2008 eResearch, ITEE, 
#                 The University of Queensland
######################################################
'''

import Pmw, sys, urllib2, urllib, string, webbrowser
from modules.idlelib.TreeWidget import TreeItem, TreeNode
from Tkinter import PhotoImage
from pymol import cmd
from urllib2 import URLError, HTTPError 
from xml.dom.minidom import parseString
from datetime import datetime
from user import home
import platform

def __init__(self):
    self.annotationService = None
    self.menuBar.addmenuitem('Plugin', 'command',
        'AnnoCryst',label = 'AnnoCryst', 
        command = lambda s=self: createAnnotationService(s))
    cmd.extend('annotate', lambda s=self : annotateFromCmd(s))
    cmd.extend('annotations', 
               lambda model, s=self : showAllAnnotations(model, s))
    cmd.extend('remoteurl', lambda pdbURL, s=self : readRemoteURL(pdbURL,s))
    cmd.extend('remotepdb', lambda pdbCode, s=self : readRemotePDB(pdbCode,s))
    
def createAnnotationService(app):
    if (app.annotationService == None):
        app.annotationService = AnnotationService(app)
    else:
        app.annotationService.dialog.show()
        
def showAllAnnotations(model, app):
    if (app.annotationService == None):
        createAnnotationService(app)
    app.annotationService.showAllAnnotations(model)
    
def annotateFromCmd(app):
    # selection, type, description, 
    if (app.annotationService == None):
        createAnnotationService(app)
    app.annotationService.annotate()

# Load a remote file: specify the full URL    
def readRemoteURL(pdbURL, app):
    if (app.annotationService == None):
        createAnnotationService(app)
    app.annotationService.openRemote(pdbURL)

# Load a remote file: specify the pdb code only
def readRemotePDB(pdbCode, app):    
    if (app.annotationService == None):
        createAnnotationService(app)
    app.annotationService.openRemoteByPDBCode(pdbCode)

class AnnotationService:
    def __init__(self,app):
        print("\nWrite the following in the PyMOL command window:")
        print("remotepdb 3ait")
        parent = app.root
        self.parent = parent
        filepath = home
        if platform.system() == 'Windows':
            if filepath == '' or filepath == None:
                filepath = "C:\\"
            if filepath[-1] != "\\":
                filepath += "\\"
        if platform.system() == 'Linux':
            if filepath == '' or filepath == None:
                filepath = "~/"
            if filepath[-1] != "/":
                filepath += "/"
        self.settingsFile = "%sannocryst-settings.xml" % filepath
        self.createSettingsDialog()
        self.createMainWindow()
        self.selection = "sele"
        if (len(cmd.get_names('selections')) == 0):
            cmd.select("all") 
        self.loadedOntology = None
        self.loadedModels = {}
        self.selectedText= ""
        self.selectedAnno = None
        self.annotationsLoaded = False

    def createMainWindow(self):
        # constructs the main AnnoCryst window
        self.dialog = Pmw.Dialog(self.parent,
            buttons = ('AnnoCryst Settings', 'Exit'),
            title = 'AnnoCryst',
            command = self.handleMainWindowButtons)
        self.notebook = Pmw.NoteBook(self.dialog.interior(), 
            raisecommand = self.refreshAnnotationView)
        self.notebook.component('hull').configure(height=380, width=400)
        self.notebook.pack(fill='both',expand=1,padx=10,pady=10)
        self.status = Pmw.ScrolledText(self.dialog.interior(),
                labelpos = 'w',
                label_text='Status: ', 
                usehullsize = 1,
                hull_width = 380,
                hull_height = 20)
        self.status.configure(text_state = 'disabled')
#        self.status.component('text').configure(relief='flat', 
#                                                background='SystemMenu')
        self.status.component('text').configure(relief='flat', 
                                                background='white')
        self.status.pack(padx = 5, pady = 5, fill = 'both', expand = 1)
        
        page = self.notebook.add('Open Model')
        self.remoteURL = Pmw.EntryField(page,
            labelpos='nw',
            label_text='URL of the model to open in AnnoCryst:',
            value = '')
        self.remoteURL.pack(fill='x', padx=4, pady=1)
        openButtonBox = Pmw.ButtonBox(page, hull_width=100, hull_height=20)
        openButtonBox.pack()
        openButtonBox.add('Open', command = self.openRemote)
        #openButtonBox.setdefault('Open')
        page = self.notebook.add("Browse Annotations")
        self.current = Pmw.ScrolledText(page, 
                usehullsize = 1,
                hull_width = 380,
                hull_height = 30)
        self.current.configure(text_state = 'disabled')
#        self.current.component('text').configure(relief='flat', 
#                                                 background='SystemMenu')
        self.current.component('text').configure(relief='flat', 
                                                 background='white')
        self.current.pack(padx = 5, pady = 5, fill = 'both', expand = 1)
        self.tree_item = AnnotationTreeItem("",isTopLevel=True)
        self.sc = Pmw.ScrolledCanvas(page,
            borderframe = 1,
            usehullsize = 1,
            hull_width = 400,
            hull_height = 270)
        self.sc.pack(fill='x',padx=4,pady=1)
        self.sc.interior().config(bg="white")
        self.node = AnnotationTreeNode(self.sc.component('canvas'), 
                                       None, self.tree_item)
        self.node.update()
        self.node.expand()
        self.browseButtonBox = Pmw.ButtonBox(page, hull_width=200, hull_height=20)
        self.browseButtonBox.add('Refresh annotations', command = self.refreshAnnotationView)
        self.browseButtonBox.add('Copy text to clipboard', command = self.copyText)
        self.browseButtonBox.add('Delete annotation', command = self.deleteAnnotation)
        self.browseButtonBox.pack()
        page = self.notebook.add('Annotate')
        group = Pmw.Group(page,tag_text='Title and Type')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5) 
        self.title = Pmw.EntryField(group.interior(),
            labelpos='w',
            label_text='Title:  ')
        self.title.pack(fill='x',padx=4,pady=1)
        self.type = Pmw.ComboBox (group.interior(),
            labelpos = 'w',
            label_text = 'Type:',
            scrolledlist_items = ('Comment', 'Rating', 'Question', \
                                  'SeeAlso', 'Feedback', 'Reference', 'Keyword'),
            selectioncommand = self.updateDescriptionUI,
            entryfield_entry_state="readonly")
        self.type.pack(fill='x',padx=4,pady=1)
        self.type.selectitem('Comment')
        self.descgroup = Pmw.Group(page,tag_text='Description')
        self.descgroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        #TODO: support for policies
        #group = Pmw.Group(page,tag_text='Access Control (optional)')
        #group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.description = Pmw.ScrolledText(self.descgroup.interior(),
            usehullsize = 1,
            hull_width = 400,
            hull_height = 100)
        self.kwdescription = Pmw.ScrolledText(self.descgroup.interior(),
            usehullsize = 1,
            hull_width = 400,
            hull_height = 30)
        self.ontologyBrowser = Pmw.ScrolledCanvas(self.descgroup.interior(),
            usehullsize = 1,
            hull_width = 250,
            hull_height = 100,
            borderframe = 1)
        self.ontologyBrowser.interior().config(bg="white")
        self.addKeywordButtonBox = Pmw.ButtonBox(self.descgroup.interior(),
		hull_width=50, hull_height=20)
        self.addKeywordButtonBox.add('Add Keyword', command = self.addKeyword)

        self.keyword = Pmw.ScrolledText(self.descgroup.interior(),
            labelpos='w',
            label_text='Keywords:  ', 
	    usehullsize =1,
	    hull_height=40,
	    hull_width=380)
        self.seeAlsoNotebook = Pmw.NoteBook(self.descgroup.interior())
        self.seeAlsoNotebook.component('hull').configure(height=100, width=400)
        seeAlsoPage = self.seeAlsoNotebook.add('External')
        self.seeAlsoExternalURL = Pmw.EntryField(seeAlsoPage,
            labelpos='w',
            label_text='URL:  ')
        self.seeAlsoExternalURL.pack(fill='x',padx=4,pady=1)
        seeAlsoPage = self.seeAlsoNotebook.add('Local')
        self.seeAlsoLocalFile = Pmw.EntryField(seeAlsoPage,
            labelpos='w',
            label_text='File:  ',
            value='File upload not yet implemented')
        self.seeAlsoLocalFile.pack(fill='x',padx=4,pady=1)
        self.updateDescriptionUI(self.type)
        self.annotateButtonBox = Pmw.ButtonBox(page, 
           hull_height=20, hull_width=200)
        self.annotateButtonBox.pack()
        self.annotateButtonBox.add('Reset', command = self.reset)
        self.annotateButtonBox.add('Annotate', command = self.annotate)
        #self.annotateButtonBox.setdefault('Annotate')
        self.parent.focus_set()
        self.annotateButtonBox.alignbuttons()
        self.dialog.show()
        
    def createSettingsDialog(self):
        # constructs a dialog for changing the settings
        self.loadSettings()
        self.settingsDialog = Pmw.Dialog(self.parent,
            buttons = ('Cancel', 'Save'),
            defaultbutton = 'Save',
            title = 'AnnoCryst Settings',
            command = self.saveSettings,
            deactivatecommand = self.saveSettings)
        attrs = self.settings.keys()
        attrs.sort()
        for att in attrs:
            entryfield = Pmw.EntryField(self.settingsDialog.interior(),
                                        labelpos='w',
                                        label_text='%s:  ' % att,
                                        value=self.settings[att],
                                        entry_width=80)
            entryfield.pack(fill='x',padx=4,pady=1)
            setattr(self,att,entryfield)
        self.settingsDialog.withdraw()
        
    def saveSettings(self, result='Cancel'):
        if result == 'Save':
            try:
            
                settingsStr = "<annocryst>\n"
                for k in self.settings.keys():
                    newvalue = getattr(self, k).getvalue()
                    self.settings[k] = newvalue
                    settingsStr += "<%s>%s</%s>\n" % (k, newvalue, k)
                settingsStr += "</annocryst>\n"
                settingsfile = open(self.settingsFile,'w')
                settingsfile.write(settingsStr)
                settingsfile.close()
                print "Settings saved in %s" % self.settingsFile
            except:
                print "Unable to save settings"
        else:
            for k in self.settings.keys():
                entryfield = getattr(self, k)
                entryfield.setvalue(self.settings[k])
        self.settingsDialog.withdraw() 
        
    def loadSettings(self):
        # default settings
        self.settings = {\
            'keywordOntologyURL': "http://maenad.itee.uq.edu.au/agerber/po.owl",
            'keywordOntologyNamespace': "http://www.proteinontology.info/po.owl",
            'annotationServerURL': "http://maenad.itee.uq.edu.au:8080/Annotea/AnnoteaServlet",
            'uploadServerURL': "http://maenad.itee.uq.edu.au:8080/Annotea/FileUploadServlet",
            'username': "Anonymous",
            'pdbRepositoryURL': "http://maenad.itee.uq.edu.au:8080/harvanapdb/au.edu.uq.itee.eresearch.harvana.gwt.Main/pdb/"
        }
        try:
            settings = open(self.settingsFile,"r")
            dom = parseString(settings.read())
            for k in self.settings.keys():
                elems = dom.getElementsByTagName(k)
                if len(elems) > 0 and len(elems[0].childNodes) > 0:
                    self.settings[k] = elems[0].childNodes[0].nodeValue
        except:
            print "Unable to read settings from %s, using AnnoCryst defaults" % self.settingsFile
            
    def handleMainWindowButtons(self,result):
        # hide or show UI dialogs
        if result=='AnnoCryst Settings':
            self.settingsDialog.show()
        else:
            self.dialog.withdraw()
            
    def openRemoteByPDBCode(self, pdbCode=''):
        pdbURL = self.pdbRepositoryURL.getvalue() + pdbCode + ".pdb"
        self.openRemote(pdbURL)
        
    # button actions for open page
    def openRemote(self, pdbURL=''):
        # load a model from a URL
        self.status.setvalue("Loading model, please wait...")
        try:
            if pdbURL == '':
                pdbURL = self.remoteURL.getvalue()
            else:
                self.remoteURL.setvalue(pdbURL)
            httpRequest = urllib2.Request(url=pdbURL)
            pdbHttpHandle = urllib2.urlopen(httpRequest) 
            pdbStr = pdbHttpHandle.read()
            modelName= string.split(string.split(str(pdbURL),"/")[-1],".")[0]
            cmd.read_pdbstr(pdbStr, modelName)
            self.status.setvalue(modelName + " loaded")
            self.loadedModels[modelName] = pdbURL
            cmd.disable('all')
            cmd.enable(modelName)
        except:
            self.status.setvalue("Unable to load model " + pdbURL)
            print "Unable to load model %s:" % pdbURL, sys.exc_info()[0]
    
    # button actions for browse page
    def copyText(self):
        copyButton = self.browseButtonBox.button("Copy text to clipboard")
        copyButton.clipboard_clear()
        copyButton.clipboard_append(self.selectedText, type='STRING')

    
    def deleteAnnotation(self):
        anno = self.selectedAnno
        if anno != None and anno != '':
            try:
                req = RequestWithMethod(anno,method="DELETE")
                response = urllib2.urlopen(req) 
                self.showAllAnnotations()
                self.status.setvalue("Annotation deleted")
            except:
                print "Unable to delete annotation"
                self.status.setvalue("Unable to delete annotation")
                
    # button actions for annotate page   
    def reset(self):
        # clears annotation fields
        self.title.setvalue('')
        self.description.setvalue('')
        self.kwdescription.setvalue('')
        self.keyword.setvalue('')
        
    def annotate(self):
        # create and post annotation of current selection
        self.annoIDs = []
        self.annoView = cmd.get_view()
        self.annoIDs = cmd.index(self.selection)
        contextStr = ''
        viewStr = ''
        contextModel = ''
        dateStr = self.makeDateString()
        for i in self.annoIDs:
            # TODO: raise an error if multiple models are selected
            # currently annotates only the model from the atom first in the selection list
            if contextModel == '':
                contextModel = i[0]
                contextStr += "%i" % i[1]
            else:
                if i[0] == contextModel:
                    contextStr += ",%i" % i[1]
        for i in self.annoView:
            if viewStr == '':
                viewStr = "%f" % i
            else:
                viewStr += ",%f" %i
        # if there is a selection or a model visible, annotate that, if not,
        # if the URL in the open page UI has been loaded, it is the most recent model, so annotate that
        uiURL = self.remoteURL.getvalue()
        annoURI = ''
        if uiURL in self.loadedModels.values():
            annoURI = uiURL  
        if (annoURI == '' or annoURI == None) and contextModel != '':
            annoURI = self.loadedModels[contextModel]            
        if annoURI != None and annoURI != '':
            annoType = self.type.getvalue()[0]
            annoCreator = self.username.getvalue()
            annoTitle = self.title.getvalue()
            annoDesc = self.description.getvalue()
            if annoType == "Keyword":
                annoDesc = self.kwdescription.getvalue()
                keywords = self.keyword.getvalue().split(",")
                anno = self.createAnnotationXML(annoURI, type = annoType, \
                    creator = annoCreator, title = annoTitle,\
                    context = contextStr, body = annoDesc, \
                    keywords=keywords, view=viewStr, date=dateStr)
            elif annoType == "SeeAlso":
                annoExtURL =  self.seeAlsoExternalURL.getvalue()
                annoLocalFile = self.seeAlsoLocalFile.getvalue()
                if annoExtURL.find("http") == 0:
                    anno = self.createAnnotationXML(annoURI, type = annoType, \
                        creator = annoCreator, title = annoTitle, date=dateStr,\
                        context = contextStr, extRef = annoExtURL, view=viewStr)
                #elif annoLocalFile != None and annoLocalFile != '':
                #TODO: do file upload
                else:
                    self.status.setvalue("Unable to create annotation: invalid external URL")
                    return
            else:
                anno = self.createAnnotationXML(annoURI, type = annoType, \
                    creator = annoCreator, title = annoTitle, date=dateStr,\
                    context = contextStr, body = annoDesc, view=viewStr)
            
            try:
                req = urllib2.Request(self.annotationServerURL.getvalue(), anno)
                response = urllib2.urlopen(req)
            except HTTPError, e:
                if e.code == 201:
                    self.showAllAnnotations()
                    self.status.setvalue("Annotation created")
                    
    def makeDateString(self):
        date = datetime.utcnow()
        dateStr = "%s-%s-%sT%s:%s:%sZ" % (date.year,\
            self.makeDatePartString(date.month), self.makeDatePartString(date.day),\
            self.makeDatePartString(date.hour), self.makeDatePartString(date.minute),\
            self.makeDatePartString(date.second))
        return dateStr
    
    def makeDatePartString(self, part):
        partStr = "%i" % part
        if part < 10:
            partStr = "0" + partStr
        return partStr
    
    def createAnnotationXML(self, annoURI, type='Comment', context='', title='', \
            creator='', language='en', created='', date='', length='', body = '', \
            view = '', keywords='', extRef = ''):
        # construct the XML for the annotation to be sent to the Annotea server
        anno = "<?xml version=\"1.0\"?>"
        anno += "<r:RDF xmlns:r='http://www.w3.org/1999/02/22-rdf-syntax-ns#'"
        anno += " xmlns:a='http://www.w3.org/2000/10/annotation-ns#'"
        anno += " xmlns:d='http://purl.org/dc/elements/1.1/'"
        anno += " xmlns:h='http://www.w3.org/1999/xx/http#'>\n"
        anno += "<r:Description>\n"
        anno += "<r:type r:resource='http://www.w3.org/2000/10/annotation-ns#Annotation'/>\n"
        if type == 'Keyword':
            anno += "<r:type r:resource='http://metadata.net/wannotea/semantic-annotation.owl#SemanticAnnotation'/>\n" 
        else:
            anno += "<r:type r:resource='http://www.w3.org/2000/10/annotationType#%s'/>\n" % type 
        anno += "<a:annotates r:resource='%s'/>\n" % annoURI
        anno += "<a:context>"
        if view != '' and view != None:
            anno += "view:%s;" % view
        if context != '' and context != None:
            anno += "ids:%s" % context
        anno += "</a:context>\n" 
        #if context != '' and context != None and view != '' and view != None:
        #    anno += "<a:context><View xmlns='http://maenad.itee.uq.edu.au/agerber/foo'>%s</View>" % repr(view)
        #    anno += "<Context xmlns='http://maenad.itee.uq.edu.au/agerber/foo'>%s</Context></a:context>" % context
        # TODO: store view data in annotation properly
        #anno += "<Context xmlns='http://maenad.itee.uq.edu.au/agerber/pymol.rdfs'>\n"
        #anno += "<Model url='%s'/>\n" % annoURI
        #if view != '' and view != None:
        #    anno += "<CameraView>%s</CameraView>\n" % view
        #anno += "</Context>\n"
        if title != '' and title != None:
            anno += "<d:title>%s</d:title>\n" % title
        if creator != '' and creator != None: 
            anno += "<d:creator>%s</d:creator>\n" % creator
        anno += "<d:language>%s</d:language>\n" % language
        if created != '' and created != None:  
            anno += "<a:created>%s</a:created>\n" % created
        if date != '' and date != None:
            anno += "<d:date>%s</d:date>\n" % date
        if body != '' and body != None:
            anno += "<a:body>\n<r:Description><h:ContentType>text/html</h:ContentType>\n"
            anno += "<h:ContentLength>%s</h:ContentLength>" % length 
            anno += "<h:Body r:parseType='Literal'>\n"
            anno += "<html xmlns='http://www.w3.org/1999/xhtml'>\n"
            anno += "<head><title>%s</title></head>\n" % title
            anno += "<body>%s</body></html></h:Body>\n" % body
            anno += "</r:Description></a:body>\n"
        if extRef != '' and extRef != None:
            anno += "<a:body r:resource=\"%s\"/>" % extRef
        if keywords !='' and keywords != None:
            for k in keywords:
                keyword = k.strip()
                # TODO: don't add the keyword if it is not from the loaded ontology
                #if keyword in self.ontology_tree_item.keywords:
                if True:
                    anno += "<term xmlns='http://metadata.net/wannotea/semantic-annotation.owl#' "
                    anno += "r:resource='%s#%s'/>" % (self.keywordOntologyNamespace.getvalue(), keyword)
                else:
                    if keyword != '':
                        print "Warning: Invalid keyword \'%s\' not added to annotation" % keyword
                    

        anno += "</r:Description></r:RDF>"
        return anno

    def refreshAnnotationView(self, page = None):
        if page == None or (page == 'Browse Annotations' and not self.annotationsLoaded):
            self.showAllAnnotations()
    def clearAnnotations(self):
        self.tree_item = AnnotationTreeItem("",isTopLevel=True)
        self.node = AnnotationTreeNode(self.sc.component('canvas'), 
                                       None, self.tree_item)
        self.node.update()
        self.node.expand()
        
    def showAllAnnotations(self, url=''):
        # get annotations for the current model, and display in tree
        # looks for an url first as a param, then from the url ui field (if it's loaded), 
        # finally from current graphical selection
        annoIDs = cmd.index(self.selection)
        self.clearAnnotations()
        if url == '' or url == None:
            uiURL = self.remoteURL.getvalue()
            if uiURL in self.loadedModels.values():
                url = uiURL
        if url == '' and len(annoIDs) > 0:
            contextModel = annoIDs[0][0]
            url = self.loadedModels[contextModel]
        if url == '' or url == None:
            self.current.setvalue("No model loaded in AnnoCryst")
            return
        self.current.setvalue("Showing annotations for: %s" % url)
        annotea = "%s?w3c_annotates=%s" % (self.annotationServerURL.getvalue(), url)
        try:
            httpRequest = urllib2.Request(url=annotea)
            httpHandle = urllib2.urlopen(httpRequest)
            #self.printData = 0
            #self.tmpAnnotation = {}
            annoteaRdfXml = httpHandle.read()
            dom = parseString(annoteaRdfXml)
            self.tree_item = AnnotationTreeItem(dom.documentElement,isTopLevel=True)
            self.node = AnnotationTreeNode(self.sc.component('canvas'), 
                                           None, self.tree_item)
            self.node.setannotationservice(self)
            self.node.update()
            # expand all nodes in tree
            self.node.expand()
            for child in self.node.children:
                child.expand()
            self.status.setvalue("Annotations loaded")
            self.annotationsLoaded = True
        except URLError, e:
            self.status.setvalue("No annotations found")
        except HTTPError, e:
            print "Unable to load annotations"
            self.status.setvalue("Unable to load annotations")


    def updateDescriptionUI(self,result):
        # changes the UI depending on the type of annotation being created
        self.keyword.pack_forget()
        self.ontologyBrowser.pack_forget()
        self.addKeywordButtonBox.pack_forget()
        self.description.pack_forget()
        self.kwdescription.pack_forget()
        self.seeAlsoNotebook.pack_forget()
        if result == 'Keyword':
            self.kwdescription.pack(fill='x',padx=4,pady=1)
            self.ontologyBrowser.pack(fill='x', padx=4, pady=1)
            self.addKeywordButtonBox.pack()
            self.keyword.pack(fill='x',padx=4,pady=1)
            self.loadOntology()
        elif result == 'SeeAlso':
            self.seeAlsoNotebook.pack(fill='both',expand=1,padx=10,pady=10)    
        else:
            self.description.pack(fill='x',padx=4,pady=1)

    def loadOntology(self):
        # retrieves the ontology and displays classes in ontology browser
        ontologyCanvas = self.ontologyBrowser.component('canvas')
        ontURL = self.keywordOntologyURL.getvalue()
        if self.loadedOntology == ontURL:
            return
        if ontURL != '' and ontURL != None:
            try:
                ontReq = urllib2.Request(url=ontURL)
                ontHandle = urllib2.urlopen(ontReq)
                ontContent = ontHandle.read()
                ontDom = parseString(ontContent) 
                self.ontology_tree_item = OntologyTreeItem(ontDom.documentElement)
                self.ontology_tree_node = OntologyTreeNode(ontologyCanvas, 
                                                    None, self.ontology_tree_item)
                self.ontology_tree_node.setannotationservice(self)
                self.ontology_tree_node.update()
                self.ontology_tree_node.expand()
                self.loadedOntology = ontURL
            except:
                print "Unable to load ontology: %s" % ontURL
                self.status.setvalue("Unable to load ontology")
                
    def addKeyword(self):
        # copies selected keyword from ontology viewer to keyword field
        kw = self.selectedKeyword
        if kw != '' and kw != None:
            oldval = self.keyword.getvalue().strip()
            if oldval == None:
               oldval = ''
            if oldval != '':
               oldval += ','
            self.keyword.setvalue(oldval + kw)
                            
    def selectAnnotation(self,context_ids, view):
        # highlights the context of the selected annotation
        select_str = "none"
        cmd.select(self.selection, 'none')
        cmd.set_view(view)
        # only add one id at a time to selection to avoid crashing PyMOL with large selection
        for id in context_ids:
            if id != "":
                select_str = self.selection + " or id %s" % id
        
                try:
                    cmd.select(self.selection, select_str)
                    cmd.indicate(self.selection)
                    self.current.setvalue("Showing context for annotation")
                except:
                    print "Unable to select context"
                    
    def deselectAnnotation(self):
        cmd.indicate("none")
        url = self.remoteURL.getvalue()
        if url == '' or url == None:
            self.current.setvalue("No model loaded in AnnoCryst")
        else:
            self.current.setvalue("Showing annotations for %s" % url)

## tree widget classes for displaying ontology keywords

class OntologyTreeNode(TreeNode):
    def __init__(self, canvas, parent, item):
        TreeNode.__init__(self, canvas, parent, item)
        self.classicon = """
R0lGODdhEAAQAOMPAAAAAIAAAACAAICAAAAAgIAAgACAgMDAwICAgP8AAAD/AP//AAAA//8A/wD/
/////ywAAAAAEAAQAAAEI/DJSau9+IXN9+2g1VFexWVkiWoqeq5sAMfi974miIvh7GMRADs=
"""
    # overload to load the icon from a string instead of from an image file
    def geticonimage(self, name):
        if name == "none":
            return None;
        if name != "Class" :
            return TreeNode.geticonimage(self, name)
        try:
            return self.iconimages[name]
        except KeyError:
            pass
        image = PhotoImage(master=self.canvas,data=self.classicon)
        self.iconimages[name] = image
        return image
    
    def setannotationservice(self, as1):
        self.annotationservice = as1

    #overload to set annotation service for children
    def draw(self,x,y):
        result = TreeNode.draw(self,x,y)
        for child in self.children:
            child.setannotationservice(self.annotationservice)
        return result

    # overload select to send selected keyword to annotation service
    def select(self, event=None):
        TreeNode.select(self, event)
        if self.annotationservice != None:
                self.annotationservice.selectedKeyword = self.item.GetText()

	
class OntologyTreeItem(TreeItem):
    def __init__(self, node, class_dict=None):
        self.children = []
    	self.node = node
    	if node.nodeType == node.ELEMENT_NODE:
    	    self.tag = node.nodeName
            if self.tag == ("rdf:RDF"):
                allclasses = node.getElementsByTagNameNS(\
                       'http://www.w3.org/2002/07/owl#','Class')
                self.class_dict = {}
                self.equiv_dict = {}
                tmpchildren = {}
                # find all OWL subclass relationships and store in class_dict
                # TODO: deal with equivalent classes
                for c in allclasses:
                    cID = c.getAttributeNS('http://www.w3.org/1999/02/22-rdf-syntax-ns#','ID')
                    cSuperElems = c.getElementsByTagNameNS(\
                                    'http://www.w3.org/2000/01/rdf-schema#','subClassOf')
                    if len(cSuperElems) > 0:
                        cSuper = cSuperElems[0].getAttributeNS(\
                                    'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                                    'resource').replace("#","")
                        if cSuper == None or cSuper == '':
                            cSuperClassElem = cSuperElems[0].getElementsByTagNameNS(\
                                        'http://www.w3.org/2002/07/owl#','Class')
                            if len(cSuperClassElem) > 0:
                                cSuper = cSuperClassElem[0].getAttributeNS(\
                                        'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                                        'about').replace("#","")
                                if cSuper == None or cSuper == '':
                                    cSuper = cSuperClassElem[0].getAttributeNS(\
                                        'http://www.w3.org/1999/02/22-rdf-syntax-ns#','ID')
                        if cSuper == None or cSuper == '':
                            # can't get the super class - might be a restriction etc
                            if cID:
                                tmpchildren[cID] = c
                        else:    
                            if not cSuper in self.class_dict:
                                self.class_dict[cSuper] = [c]
                            else:
                                self.class_dict[cSuper].append(c)
                    else: 
                        if cID:
                            tmpchildren[cID] = c
                # save list of keywords
                self.keywords = tmpchildren.keys()
                # iterate over all classes that have super classes, 
                # remove them from the children of the top element
                for childList in self.class_dict.values():
                    for child in childList:
                        childId = child.getAttributeNS(\
                            'http://www.w3.org/1999/02/22-rdf-syntax-ns#','ID')
                        if childId == None or childId == '':
                            childId = child.getAttributeNS(\
                                'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                                'about').replace("#","")
                        if childId in tmpchildren:
                            del tmpchildren[childId]
                self.children = tmpchildren.values()
            elif self.tag.find("Class") != -1:
                # look up the subclasses (children) of this class
                self.class_dict = class_dict
                cID = self.GetText()
                if cID in class_dict:
                    self.children = class_dict.get(cID)
                
    def GetLabelText(self):
        if self.tag == "rdf:RDF":
            return "Select keywords:"
        else:
            return ""
    
    def GetText(self):
        # return the ID of the class that this node represents
        if self.tag == "rdf:RDF":
            return " "
        text = self.node.getAttributeNS('http://www.w3.org/1999/02/22-rdf-syntax-ns#','ID')
        if text == None or text == '':
            text = self.node.getAttributeNS(\
                'http://www.w3.org/1999/02/22-rdf-syntax-ns#','about').replace("#","")
        return text

    def IsExpandable(self):
        return len(self.children) > 0
        
    def GetSubList(self):
        parent = self.node
        itemlist = [OntologyTreeItem(node,self.class_dict) for node in self.children]
        return itemlist   
    
    def GetIconName(self):
        if self.tag == "rdf:RDF":
            return "none"
        else:
            return "Class"

## tree widget classes for displaying annotations           
class AnnotationTreeNode(TreeNode):
    
    def __init__(self, canvas, parent, item):
        TreeNode.__init__(self, canvas, parent, item)
        self.annotationservice = None;
        self.comment_icon = """
R0lGODdhEAAQAOeZAERERFJSUkxaZ1hbXU1eb1xcXF1dXVJgbmVlZVRthWtra3BwcG14g2l+k3Z7
gHN9h3V9hGGFq2uFnnaMom6QsXSQq5CQkJGRkXyZtXObw4GZsXyevpmZmYmiuYakw5+fn4elwoqk
wYWmxJCluKGhoYmpxpmntZynspKvz5K33LGxsbW1tam4yqm5yKC71p2/3qq/07u7u6vA1qDD5J/E
6qnD3KXF5LDG3rnFz7LH3bjG1LXH2a3K5cXFxb/J0qrO86/N7bfN48nT3czX4rjd/73c/Lre/9fZ
3MTh+7/j/8Xj/8Dn/97e3sbl/9Hi9srk/8fo/83m/+Di49zj6svo/83n/+Li4s7o/9rl8OPj483q
/8/p/9Dq/9Pp/+Xl5dHr/9Xq/9Tr/+Do8OLo7dTs/9br/9nt/93s++jr7tzu/+vr69rw/97u/+js
7+zr6+zs7Ojt8t3x/+Hw+uzt79/x/+ru8t7y/+jv9unv9Orv9Obx/O/v7+zw9ePz/+Xy//Dw8O3x
9ef0//Ly8uf2/+T4//Pz8+f4/+v3/+n5/+33//b19ff39/n5+Pr59vn5+e/+//j6/Pr6+vv7+/b/
//78+vj///z+//v///3/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////ywAAAAAEAAQAAAI3wAzCYRU
p42iRozyjEETSaBDOHMoOXTIR8gbh2LcTNwISAeTTHekWLq00SEmLCasDMEzKVHJTJUQyZHBwYee
Q34CDbLzKFAgQ4bMrMnhAMeZPmzS0FlKxxChLlyu3IDAwkmcMmDCaNVa5QkVJC4GwNhB5kuULVDS
QlGipAkPDACmjAiihYpatUmAgDiAQFKLDjWKGBlshIgRGx4SGFCRSdCJCiFQvPixZIYICgQKfLiY
6U+MBw0mlKCxQUCABSrUTHTkpQcJDSkkKFhhZdHLSIWyZIhwQfXLiUcYWNgzMSAAOw==
"""
        self.question_icon = """
R0lGODdhEAAQAOfAAAAASQAAVwAAXAAAZAAfcAAdjQ80iBE6nQ87pyE+hRlGpURERDZOkFJSUkxa
Z1hbXU1eb1xcXF1dXVJgbmVlZVRthWtra01xuXBwcEd1yG14g156q2l+k3Z7gHN9h3V9hGGFq2uF
nmCJy2qIvXaMom6QsXSQq5CQkJGRkXyZtXObw4GZsX+ZxYWYvXyevoSaxZmZmYKgzomiuYakw5+f
n4elwoqkwYWmxJCluKGhoYmpxoip05mntZynspenwpKvz52z0pK33JO247GxsZm23LW1tam4yqm5
yKC71p2/3p+/6Kq/07u7u6vA1qDD5J/E6qnD3KXF5LDG3rLF4rnFz7LH3bjG1LXH2a3K5cXFxb/J
0qrO86/N7bfN47bR8bjS7b7T68nT3czX4rjd/73c/Lre/9fZ3Lzh/8Th+7/j/8Xj/8Dn/97e3sbl
/9Hi9srk/8fo/83m/+Di49zj6svo/83n/+Li4s7o/9rl8OPj483q/8/p/9Dq/8/r/9Pp/+Xl5dHr
/9Xq/9Tr/+Do8OLo7c7v/9Ts/9br/9Lv/9nt/93s+9bv/+jr7tzu/+vr69rw/97u/+js7+zr6+zs
7Ojt8tf0/93x/+Hw+uzt79/x/+ru8t7y/+jv9unv9Orv9Obx/O/v7+zw9ePz/+Xy//Dw8O3x9ef0
//Ly8uf2/+T4//Pz8+f4/+v3/+n5/+33/+b7//b19ff39+n///n5+Pr59vn5+e/+//j6/Pr6+vv7
+/b///78+vj///z+//v///3/////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////ywAAAAAEAAQAAAI+wCBCbyl
KRIsWrM8EWKES6BDSph2OXQYKswkh4MkOXxxgYEPgaWssAHGSU4vX8BaCJiSIAAQYL/w8LAjppMu
V8BYFACzYcAIXq0uNYGh5ROrUaZQbbJlAMGXVYkeVelARZEoSI0y8YqhoE8qP3zuSPlgxI2lQ4EE
vRKhRFadN3TQIHmw5IohQHH2LCLiBZEaNW2wpFgwB0cXPXQq7QBA4AycNFxqTKCQ64gMKGQKCTmQ
YUyZKDMqSBgC7FQPEzZ+JNmyxsmNEhAi0LgIjBQTDxxI6HjiwkEDDEMcTaz1J0uOFUFCWChiJ9ZE
h7hU5VEBAoXw59jNaDgBamJAADs=
"""
        self.semantic_icon = """
R0lGODdhEAAQAOeHAAoYbERERFJSUkxaZ1hbXU1eb1xcXF1dXVJgbmVlZVRthWtra3BwcG14g2l+
k3Z7gHN9h3V9hGGFq2uFnnaMom6QsXSQq5CQkJGRkXyZtXObw4GZsXyevpmZmYmiuYakw5+fn4el
woqkwYWmxJCluKGhoYmpxpmntZynspKvz5K33LGxsbW1tam4yqm5yKC71p2/3qq/07u7u6vA1qDD
5J/E6qnD3KXF5LDG3rnFz7LH3bjG1LXH2cXFxb/J0qrO86/N7bfN48nT3czX4rjd/73c/Lre/9fZ
3L/j/8Xj/8Dn/97e3sbl/9Hi9srk/8fo/+Di49zj6svo/+Li4trl8OPj483q/8/p/9Dq/9Pp/+Xl
5dTr/+Do8OLo7dbr/93s++jr7uvr6+js7+zr6+zs7Ojt8t3x/+zt79/x/+ru8ujv9unv9Orv9Obx
/O/v7+zw9eXy//Dw8O3x9ef0//Ly8uf2//Pz8+f4/+n5//b19ff39/n5+Pr59vn5+e/+//j6/Pr6
+vv7+/78+vj///z+//v///3/////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////ywAAAAAEAAQAAAI3AAPCfyT
RkwePnvYdAEDSKDDMmcEOXT4RggZh1zGTNwoZ8eSQ2qgECq00aEhKiemDFkDoOWhlgBeDsIDYEYH
H20AwJlTB4CfOXNgAtDxIMcXmGhaorkDIAsWADgitGhixguALVe3YHUiBcALAjF4wLwC8wmAJEwA
ZAgQhUQQK1KeyJ2LBEgIBAkCufBgo4iRv0aIGLnxQcGBFYfooLAgIgWMH0pojKhQwACIi4fiyIDg
gIKJGhwGCGCwIszEPlp6lNigYsICFlP0lDwEyE4VDRIwmJ498UiDC24mBgQAOw==
"""
        self.seealso_icon = """
R0lGODdhEAAQAOezABAoQBAoUCs6WCk+YTBIYDlIZkRSbUJVZ0BYYEVTbklXZkRZbUlXdTBocEBo
UDB4QEBocE5hdFJfe1JgeyCIMFNhfFNhfVBogFdndVBwYFBokFB4YCCYQCCQgGRxfkCAkCCgQGR0
gzCYUEB40ECIgGd1giCgYGB4kGd4iDCQoDCgUG13i2CAgECYYECgQFCA0G9/jlCYYGmCm3CAoGuE
nVCI4ECgkGCQkGCQoECY0ECwUGCI0GCYgECwYGCI4HCQoHyJpX6LmX2KpmCQ4ISMnmCwQHCQwHCg
cGCY0IaQmWCgsICQsHCQ4IeQomCY8ICYoGC4UGCY/3CogHCgsICYsHCY4I6VpXCY8I+Yn4yXr4Cg
oIqYs4uYtIacsmDIMIyatZGbpWDIQIqfs4qftJCctZado3DAUJScsnCg/5SeqJWcspKet5efpZWg
q3DIQJSjspajsJ2ippWksoC4kJCowJ+nrpSqwZ6ovKCowJ+qvpC4sKGsw5+vvqqssqCvv6CwwKCw
0IDYYKmzvaq0vqq1vqu1v4DYgLC1ubC4v5DYcKy40ZDQoLK5v7C40LW6vpDQsJDgcJDYoLa+xbDA
0JDggLq/xLy/wr7Bw8HBw6DggL/EybDI4KDgkMHGy8DHy8jJxqDosMDQ4MDQ8M/PztDS09DY4LDw
wNDY8NjZ2Nra2Nzc2t/c19Dg8NDw0ODo/9D44ND48PDw//D4////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////ywAAAAAEAAQAAAI+wBnzUql
SdCbM3ccjRLIUCCpQXoeRaJk5gYYSw1nqSI0x5CpV1ISFfmQ5FPDTlq8QDLVChYoHjlIxGlYCAeU
QJwyLVLipIYIGA3lMIkRxs2RKUiGpADBoKEaNCxc6LhRZUgNExwUNPRz5UqGFj1s+Hjh4EGIhpWM
RBmyQYWKDiMgIMDSEFWaHUae/DhBoYEGIWO2wDkkkBGeUrJiTZpxgQ5ixUD6zOLDKogAAlQA0Wlk
50CENn9WzOoyawAmRbNkuTolCsUlMaEKzMoTKguNNQEA6F4gg8weKhhmefoCKNSmJRWscCl+XEIZ
gZKaWEhQos4qREQmGPDAZlZAADs=
"""
    def setannotationservice(self, as1):
        self.annotationservice = as1
        
    # overload to load the annotation icon data from strings
    def geticonimage(self, name):
        if name == "none":
            return None;
        if name != "Comment" and name != "Question" and name != "Semantic" \
        and name != "SeeAlso" and name != "Reference" and \
        name != "Feedback" and name != "Rating":
            return TreeNode.geticonimage(self, name)
        try:
            return self.iconimages[name]
        except KeyError:
            pass

        icondata = self.comment_icon
        if name == "Question":
            icondata = self.question_icon
        elif name == "Semantic":
            icondata = self.semantic_icon
        elif name == "SeeAlso":
            icondata = self.seealso_icon
        image = PhotoImage(master=self.canvas,data=icondata)
        self.iconimages[name] = image
        return image

    def select(self, event=None):
        TreeNode.select(self, event)
        text = self.item.GetText()
        if self.annotationservice != None:
            # notify the selected annotation id - for delete function
            self.annotationservice.selectedAnno = self.item.annoID
            # notify annotationservice of text - for copy text function
            self.annotationservice.selectedText = text
            # notify annotationservice to highlight context in graphical view
            if self.item.label == "context":
                contextsplit = string.split(text,";")
                if len(contextsplit) > 0:
                    view = contextsplit[0]
                    view = string.replace(view,"view:","")
                    view = string.split(view,",")
                    viewfloats = []
                    for v in view:
                        viewfloats.append(float(v))
                if len(contextsplit) > 1:
                    context = contextsplit[1]
                    context = context.replace("ids:","")
                    contextids = string.split(context,",")
                #print "context: " + repr(contextids)
                #print "view: " + repr(view)
                self.annotationservice.selectAnnotation(contextids, view)
            else:
                self.annotationservice.deselectAnnotation()
        if self.item.label == "body":
            # launch URL of SeeAlso body in browser
            if text.find("http") == 0:
                webbrowser.open(text)

    
    #overload to set annotation service for context highlighting
    def draw(self,x,y):
        result = TreeNode.draw(self,x,y)
        for child in self.children:
            child.setannotationservice(self.annotationservice)
        return result
        
class AnnotationTreeItem(TreeItem):
    def __init__(self, annotation, isTopLevel=False, label=None, id=None):
        if annotation != "":
            self.anno = annotation
        else:
            self.anno = ""
        self.isTopLevel = isTopLevel
        self.label = label
        self.isLeaf = False
        if label == None and not isTopLevel:
            self.label = self.GetText()
        if self.label == "context" or self.label=="date" or self.label == "creator" \
                or self.label == "created" or self.label == "title" \
                or self.label == "identifier" or self.label == "language":
            self.isLeaf = True
        self.annoID = id
        if annotation != "" and annotation.nodeType == annotation.ELEMENT_NODE and \
        annotation.nodeName == "rdf:Description":
            self.annoID = annotation.getAttributeNS('http://www.w3.org/1999/02/22-rdf-syntax-ns#','about')
                      
            
    def GetText(self):
        node = self.anno
        if node == "":
            return " "
        elif self.isLeaf:
            text = "" 
            for child in node.childNodes:
                if child.nodeType == node.TEXT_NODE:
                    text += child.nodeValue
            return text
        elif self.label == "body" or self.label == "term":
            bodyurl = self.anno.getAttributeNS('http://www.w3.org/1999/02/22-rdf-syntax-ns#','resource')
            # get body content if it's a resource stored on the Annotation Server, otherwise show url
            if bodyurl.find('AnnoteaServlet') != -1:
                return self.getAnnotationBody(self.anno)
            else:
                return bodyurl
        else:
            if node.nodeType == node.ELEMENT_NODE:
                nName = node.nodeName
                if nName == "rdf:RDF":
                    nName = "Annotations"
                elif nName == "rdf:Description":
                    annoType = self.GetIconName()
                    if annoType != 'none':
                        nName = self.GetIconName() + " Annotation"
                    else:
                        nName = "Annotation"
                return nName
            elif node.nodeType == node.TEXT_NODE:
                return node.nodeValue

        
    def IsExpandable(self):
        if self.isLeaf or self.label == "body":
            return False
        if self.anno == "":
            return self.isTopLevel == 1
        else:
            return self.anno.hasChildNodes()
        
    def GetSubList(self):
        if self.anno == "" or self.isLeaf:
            return None
        parent = self.anno
        children = parent.childNodes
        prelist = [AnnotationTreeItem(node,id=self.annoID) for node in children \
                   if node.nodeName != "rdf:type" and node.nodeName != "policy" \
                   and node.nodeName != "annotates" and node.nodeName != "language"]
        itemlist = [item for item in prelist if item.GetText().strip()]
        return itemlist

    def GetIconName(self):
        if (self.label == None) or (self.label == "term") or \
        (self.label == "body") or (self.label == "title") or \
        (self.label == "creator") or (self.label == "context") or \
        (self.label == "created") or (self.label == "description") or\
        (self.label == "identifier") or (self.label == "date"):
            return "none"
        node = self.anno        
        for typeNode in node.getElementsByTagNameNS(\
                            'http://www.w3.org/1999/02/22-rdf-syntax-ns#','type'):
            typeStr = typeNode.getAttributeNS(\
                'http://www.w3.org/1999/02/22-rdf-syntax-ns#','resource').replace(\
                "http://www.w3.org/2000/10/annotationType#","")
            if typeStr == "Question" or typeStr == "Rating" or \
            typeStr == "SeeAlso" or typeStr == "Feedback" or typeStr=="Reference":
                return typeStr
            if typeStr == "http://metadata.net/wannotea/semantic-annotation.owl#SemanticAnnotation":
                return "Semantic"
        return "Comment"
    
    def GetLabelText(self):
        if self.label != "Annotation":
            return self.label

    def getAnnotationBody(self, node):
        bodyContentStr = ""
        try:
            bodyReq = urllib2.Request(url=node.getAttributeNS(\
                'http://www.w3.org/1999/02/22-rdf-syntax-ns#','resource'))
            bodyHandle = urllib2.urlopen(bodyReq)    
            bodyContent = bodyHandle.read()
            bodyDom = parseString(bodyContent)
            bodyElem = bodyDom.getElementsByTagNameNS("http://www.w3.org/1999/xhtml","body")
            for body in bodyElem:
                for child in body.childNodes:
                    if child.nodeType == node.TEXT_NODE:
                        bodyContentStr += child.nodeValue
        except:
            print "Unable to read annotation body"
        return bodyContentStr.strip()
        
# override urllib2 Request to support HTTP DELETE request
class RequestWithMethod(urllib2.Request):

    def __init__(self, url, data=None, headers={}, origin_req_host=None, 
                 unverifiable=False, method=None):
        urllib2.Request.__init__(self, url, data, headers, origin_req_host, unverifiable)
        self.method = method

    def get_method(self):
        if self.method == None:
            if self.data != None:
                return "POST"
            else:
                return "GET"
        else:
            return self.method 
