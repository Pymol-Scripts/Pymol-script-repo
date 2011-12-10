## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#########################################################################
#
# Date: Oct. 2002  Author: Daniel Stoffler
#
# Copyright: Daniel Stoffler and TSRI
#
#########################################################################


import string, re, urllib, types


class CGIForm:
    """this object stores the forms data parsed from a HTML file by the
    ParseForForms class."""

    def __init__(self, url=None, name=None, method=None, action=None,
                 enctype=None, input=None, radiobutton=None, checkbutton=None,
                 select=None, textarea=None, arguments=None, hiddenInput=None,
                 fieldOrder=None):
        
        # note: something i don't understand: in the constructor, i cannot say
        # radiobutton={} and then self.radiobutton = radiobutton
        # because every new instance of this class would point to the
        # dictionary of the first instance (tested with python1.5 to 2.2).
        # Instead, i have to initialise radiobutton=None and do the test below
        
        self.url = url       # url: e.g. 'http://www.google.com'
        self.name = name     # name of this form object

        if method is None: method = 'get'
        assert string.lower(method) in ['post', 'get']
        self.method = method # default cgi method, if not specified  

        self.action = action   # the cgi action. In case of google: '/search'
        self.enctype = enctype # enctype of CGI Form

        if input is None: input = {}
        self.input = input # dictionary with all the 'inputs' of this form
                           # note that input of type 'hidden' is not added
                           # to this dict, but to self.hiddenInput, input of
                           # type 'radio' and 'checkbox' are added to
                           # different dicts too as well as input of type
                           # 'select' and 'textarea'

        if radiobutton is None: radiobutton = {}
        self.radiobutton = radiobutton  # dict holding the radiobuttons

        if checkbutton is None: checkbutton = {}
        self.checkbutton = checkbutton  # dict holding the checkbox buttons

        if select is None: select = {}
        self.select = select # dictionary with all the 'selects' of this form  

        if textarea is None: textarea = {}
        self.textarea = textarea # dictionary with all the 'textareas'

        if arguments is None: arguments = {}
        self.arguments = arguments # dict with the arguments that are used to
                                   # generate the string to be sent to the
                                   # server

        if hiddenInput is None: hiddenInput = {}
        self.hiddenInput = hiddenInput # dictionary with all the inputs of
                                       # type 'hidden'

        if fieldOrder is None: fieldOrder = []
        self.fieldOrder = fieldOrder # stores tuples of (inputname, dictname)
                                     # in the order they appear in the <FORM>


    def run(self):
        args = self.arguments.copy()
        args.update(self.hiddenInput)
        args = urllib.urlencode(args)

        if self.method == 'post':
            f = urllib.urlopen(self.url+self.action, args)
        else: # method is 'get'
            if self.action is not None:
                link = self.url+self.action+'?'+args
            else:
                link = self.url
            f=urllib.urlopen(link)
        return f.read()


    def setArguments(self, **kw):
        for k,v in kw.items():
            if type(v) in (types.FloatType, types.IntType,
                           types.LongType, types.StringType):
                self.arguments[k] = v
            else:
                pat = re.compile("[\[,\]\012]")
                arrayStr = re.sub(pat, '', str(Numeric.array(v).ravel()) )
                c = string.split(arrayStr)
                c = map( float, c )
                c = re.sub(pat, '', str(c) )
                self.arguments[k] = c


    def setURL(self, url):
        if url is None: return
        self.url = url
        if self.action is not None:
            if self.url[-1] != '/' and self.action[0] != '/':
                    self.action = '/' + self.action

    def dump(self):
        print '............'+repr(self)+'............'
        for k,v in self.__dict__.items():
            print '%s\t:'%(k,), v
        print '...............................................................'


    def getCreationSourceCode(self):
        url = self.url
        if url is None: url = ''
        name = self.name
        if name is None: name = ''
        method = self.method
        if method is None: method = 'get'
        enctype = self.enctype
        if enctype is None: enctype = ''
        action = self.action
        if action is None: action = ''
        input = self.input
        if input is None: input = {}
        radiobutton = self.radiobutton
        if radiobutton is None: radiobutton = {}
        checkbutton = self.checkbutton
        if checkbutton is None: checkbutton = {}
        select = self.select
        if select is None: select = {}
        textarea = self.textarea
        if textarea is None: textarea = {}
        arguments = self.arguments
        if arguments is None: arguments = {}
        hiddenInput = self.hiddenInput
        if hiddenInput is None: hiddenInput = {}
        fieldOrder = self.fieldOrder
        if fieldOrder is None: fieldOrder = []
        
        txt = 'CGIForm(url="'+url+'", name="'+name+'",\n'+\
              'method="'+method+'", enctype="'+enctype+'", action="'+\
              action+'",\n'+\
              'input='+`input`+',\n'+\
              'radiobutton='+`radiobutton`+',\n'+\
              'checkbutton='+`checkbutton`+',\n'+\
              'select='+`select`+',\n'+\
              'textarea='+`textarea`+',\n'+\
              'arguments='+`arguments`+',\n'+\
              'hiddenInput='+`hiddenInput`+',\n'+\
              'fieldOrder='+`fieldOrder`+')\n'

        return txt



class ParseForCGIForms:
    """ parsing HTML files for forms. The forms description of forms
    found in the html page is stored in the dictionary self.data """

    def __init__(self):
        self.data = [] # stores all found forms as objects


    def getAllTags(self, html): # currently not used
        """ parses html page for tags, merges multi-line tags in a single
        line, and puts each tag in a new line, and discards everything
        between tags"""
        
        starttagPat = re.compile('<')
        endtagPat = re.compile('>')
        newhtml = []
        tagData = ""
        for line in html:
            i = 0
            n = len(line)

            while i < n:
                startmatch = starttagPat.search(line, i)
                endmatch = endtagPat.search(line, i)
                if startmatch and endmatch:
                    if startmatch.start() > endmatch.start():
                        tagData = tagData + " " + line[:endmatch.end()]
                        newhtml.append(tagData)
                        tagData = ""
                        i = endmatch.end()
                    else:
                        tagData = line[startmatch.start():endmatch.end()]
                        newhtml.append(tagData)
                        tagData = ""
                        i = endmatch.end()
                    
                elif startmatch:
                    tagData = line[startmatch.start():-1]+ " "
                    i = n
                    break

                elif endmatch:
                    tagData = tagData + " " + line[:endmatch.end()]
                    newhtml.append(tagData)
                    tagData = ""
                    i = endmatch.end()
                else:
                    i = n
                    if line is not None and line != '' and line != '\n':
                        tagData = tagData + line[:-1]
                    break

        return newhtml


    def orderHTML(self, html):
        """ parses html page for tags, merges multi-line tags in a single
        line, and puts each tag in a new line, as well as raw data"""
        _openTag = 0
        starttagPat = re.compile('<')
        endtagPat = re.compile('>')
        newhtml = []
        tagData = ""

        for line in html:
            if line[-1] == '\n': line = line[:-1] # get rif of '\n'
            i = 0
            n = len(line)

            while i < n:
                startmatch = starttagPat.search(line, i)
                endmatch = endtagPat.search(line, i)

                if startmatch and endmatch:
                    if startmatch.start() > endmatch.end():
                        _openTag = 0
                        tagData = tagData + " " + line[i:endmatch.end()]
                        if tagData != '': newhtml.append(tagData)
                        tagData = ""
                        i = endmatch.end()
                    else:
                        tagData = line[startmatch.start():endmatch.end()]
                        rawData = string.strip(line[i:startmatch.start()])
                        rawmatch = endtagPat.search(rawData)     
                        if not rawmatch and rawData != '':
                            newhtml.append(rawData)
                        i = endmatch.end()
                        if tagData != '': newhtml.append(tagData)
                    
                elif startmatch:
                    _openTag = 1
                    rawData = string.strip(line[i:startmatch.start()])
                    if rawData != '': newhtml.append(rawData)
                    tagData = line[startmatch.start():]
                    break                 

                elif endmatch:
                    _openTag = 0
                    tagData = tagData + " " + line[i:endmatch.end()]
                    if tagData != '': newhtml.append(tagData)
                    tagData = ""
                    i = endmatch.end()
                    
                else:
                    if _openTag:
                        tagData = tagData + line
                    else:
                        rawData = string.strip(line[i:])
                        if rawData != '': newhtml.append(rawData)
                    break

        return newhtml


    def get_starttag(self, html, tag):
        data = []
        starttagPat = re.compile('<'+tag, re.IGNORECASE)

        for line in html:
            st = starttagPat.search(line)
            if st:
                data.append(line)
        return data
        

    def get_startendtag(self, html, tag):
        tagdata = []
        data = []
        _tagfound = 0 
        text = ""

        starttagPat = re.compile('<'+tag, re.IGNORECASE)
        endtagPat = re.compile('</'+tag+'>', re.IGNORECASE)

        for line in html:
            st = starttagPat.search(line)
            et = endtagPat.search(line)

            if st and not _tagfound:
                _tagfound = 1
                tagdata.append(line)
            elif et and _tagfound:
                tagdata.append(line)
                data.append(tagdata)
                _tagfound = 0
                tagdata = []
            elif _tagfound:
                tagdata.append(line)
            else:
                continue

        return data


    def get_Attr(self, html, pattern, lc=0):
        attrPat =re.compile(
            ('([%s]*=[%s]*' % (string.whitespace, string.whitespace))
            + r'(\'[^\']*\'|"[^"]*"|[-a-zA-Z0-9./:+*%?!\(\)_#=~]*))?')
        
        found = pattern.search(html)
        if found:
            attrP = attrPat.search(html[found.end():])
            if attrP:
                attr = attrP.group()[1:]
                if attr and attr[0] == '"' and attr[-1] == '"':
                    if len(attr)>2:
                        attr = attr[1:-1]
                if lc:
                    return string.lower(attr)
                else:
                    return attr
        else:
            return None
        

    def parse(self, url, html):
        """ call this method to parse html for forms """
        attrPat = re.compile('".*?"')
        namePat = re.compile('name', re.IGNORECASE)
        methodPat = re.compile('method', re.IGNORECASE)
        actionPat = re.compile('action', re.IGNORECASE)
        enctypePat = re.compile('enctype', re.IGNORECASE)

        typePat = re.compile('type', re.IGNORECASE)
        valuePat = re.compile('value', re.IGNORECASE)
        checkedPat = re.compile('checked', re.IGNORECASE)
        sizePat = re.compile('size', re.IGNORECASE)
        maxlengthPat = re.compile('maxlength', re.IGNORECASE)
        buttonstatusPat = re.compile('checked', re.IGNORECASE)

        opentagPat = re.compile('<')

        rowsPat = re.compile('rows', re.IGNORECASE)
        colsPat = re.compile('cols', re.IGNORECASE)

        # first, extract all tags in this html page, and merge multi-line
        # tags in a single line, and put every tag in a single line
        ohtml = self.orderHTML(html)

        # now lets find all tags that belong to a cgi form
        formsdata = self.get_startendtag(ohtml, tag='form')

        formnameindex = 0   # these indices (see below) are used to generate
                            # names if no name was specified in the tag 
        
        
        for form in formsdata:
            formObject = CGIForm() # this class stores the information we are
                                   # going to parse

            formObject.url = url

            inputnameindex = 0  # see formnameindex
            selectnameindex = 0
            textnameindex = 0

            # first, we process what is in side <FORM....>
            formtag = self.get_starttag(form, tag='form')[0]
            method = self.get_Attr(formtag, methodPat, lc=1)
            if method:
                formObject.method = method

            action = self.get_Attr(formtag, actionPat)
            if action:
                formObject.action = action

            enctype = self.get_Attr(formtag, enctypePat, lc=1)
            if enctype:
                formObject.enctype = enctype

            formname = self.get_Attr(formtag, namePat)
            if formname is None:
                formname = 'form_' + `formnameindex`
                formnameindex = formnameindex + 1
            formObject.name = formname

            # now, let's find all the <INPUT...> tags:
            inputsdata = self.get_starttag(form, tag='input')
            if inputsdata != []:
                radiobuttonDict = {}
                checkbuttonDict = {}
                for ipt in inputsdata:
                    inputname = self.get_Attr(ipt, namePat)
                    if inputname is None:
                        inputname = 'input_'+`inputnameindex`
                        inputnameindex = inputnameindex + 1

                    value = self.get_Attr(ipt, valuePat)
                    type = self.get_Attr(ipt, typePat, lc=1)


                    if type == 'hidden':
                        dict = {}
                        dict['value'] = value
                        formObject.hiddenInput[inputname] = dict
                        formObject.fieldOrder.append((inputname,
                                                      'hiddenInput'))
                        continue

                    if type == 'radio':
                        check = buttonstatusPat.search(ipt)
                        if check:
                            entry = (value,'on')
                        else:
                            entry = (value,'off')
                        
                        if not radiobuttonDict.has_key(inputname):
                            radiobuttonDict[inputname] = []
                            formObject.fieldOrder.append((inputname,
                                                      'radiobutton'))
                        radiobuttonDict[inputname].append(entry)
                        continue

                    if type == 'checkbox':
                        check = buttonstatusPat.search(ipt)
                        if check:
                            entry = (value, 'on')
                        else:
                            entry = (value, 'off')
                        
                        if not checkbuttonDict.has_key(inputname):
                            checkbuttonDict[inputname] = []
                            formObject.fieldOrder.append((inputname,
                                                          'checkbutton'))
                        checkbuttonDict[inputname].append(entry)
                        continue

                    inputDict = {}
                    inputDict['type'] = None
                    inputDict['value'] = None
                    inputDict['checked'] = None
                    inputDict['size'] = None
                    inputDict['maxlength'] = None

                    if type:
                        inputDict['type'] = type

                    if value:
                        inputDict['value'] = value

                    checked = self.get_Attr(ipt, checkedPat)
                    if checked:
                        inputDict['checked'] = checked

                    size = self.get_Attr(ipt, sizePat)
                    if size:
                        inputDict['size'] = size

                    maxlength = self.get_Attr(ipt, maxlengthPat)
                    if maxlength:
                        inputDict['maxlength'] = maxlength

                    formObject.input[inputname] = inputDict
                    formObject.fieldOrder.append((inputname,'input'))

                if radiobuttonDict != {}:
                    formObject.radiobutton = radiobuttonDict
                
                if checkbuttonDict != {}:
                    formObject.checkbutton = checkbuttonDict

            # now, let's find if we got <SELECT> </SELECT> tags: 
            selectdata = self.get_startendtag(form, tag='select')
            if selectdata != []:
                for select in selectdata:
                    selectDict = {}
                    options = []
                    selname = self.get_Attr(select[0], namePat)
                    if selname is None:
                        selname = 'selection_'+`selectnameindex`
                        selectnameindex = selectnameindex + 1

                    for sel in select:
                        opt = opentagPat.search(sel)
                        if not opt:
                            options.append(sel)
                        
                        size = self.get_Attr(sel, sizePat)
                        if size:
                            selectDict['size'] = size

                    selectDict['options'] = options
                                     
                    formObject.select[selname] = selectDict
                    formObject.fieldOrder.append((selname,'select'))

            # and now let's see if we got any <TEXTAREA> </TEXTAREA> tags
            textareadata = self.get_startendtag(form, tag='textarea')
            
            if textareadata != []:
                for textarea in textareadata:
                    textareaDict = {}
                    areatext = []
                    textname = self.get_Attr(textarea[0], namePat)
                    if textname is None:
                        textname = 'textarea_' + `textnameindex`
                        textnameindex = textnameindex + 1

                    for area in textarea:
                        txt = opentagPat.search(area)
                        if not txt:
                            areatext.append(area)
                    textareaDict['text'] = areatext

                    rows = self.get_Attr(textarea[0], rowsPat)
                    if rows:
                        textareaDict['rows'] = rows

                    cols = self.get_Attr(textarea[0], colsPat)
                    if cols:
                        textareaDict['cols'] = cols

                    formObject.textarea[textname] = textareaDict
                    formObject.fieldOrder.append((textname,'textarea'))

             # append this new form to self.data
            self.data.append(formObject)
            #formObject.dump()


class ParseHTML:
    """ Parses HTML files by specifying a parser.
    Current parsers: - ParseForCGIForms
    """
    
    def __init__(self, mode=None):
        self.mode = mode
        self.parsers = {} # dictionary storing the various parser objects
        # add new parsers here:
        self.parsers['forms'] = ParseForCGIForms()
                

    def parse(self, url, html):
        if html is None or html == [] or html == '':
            return

        if self.parsers.has_key(self.mode):
            self.currentParser = self.parsers[self.mode]
        else:
            return

        self.currentParser.parse(url, html)

        return self.currentParser.data


    def doit(self, url):
        """ call this with url to parse html."""
        f = urllib.urlopen(url)
        data = f.read()
        f.close()
        parsedData = self.parse(url, [data])
        return parsedData


def test1():
    f = open('../../web/example7.html','r')
    txt = f.readlines()
    P = ParseHTML(mode='forms')
    forms = P.parse(txt)
    f.close()
    return forms


def test2():
    #url = 'http://www.google.com'
    url = 'http://www.amazon.com'
    P = ParseHTML(mode='forms')
    forms = P.doit(url)
    return forms
    
if __name__=='__main__':
    f=test2()
