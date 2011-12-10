## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/InputForm/Tk/Tests/test_InputForm.py,v 1.7 2007/12/04 21:28:04 vareille Exp $
#
# $Id: test_InputForm.py,v 1.7 2007/12/04 21:28:04 vareille Exp $
#

"""
This module implements a set of functions to test the InputForm class
"""

import sys
from mglutil.regression import testplus
import numpy.oldnumeric as Numeric
import Tkinter, Pmw
from mglutil.gui.InputForm.Tk.gui import InputFormDescr, InputForm,\
     CallBackFunction
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ListChooser
from mglutil.util.misc import ensureFontCase

master = None
root = None
descr = None
Value = 0


def setUp():
    global master
    global root
    global Value
    master = Tkinter.Tk()
    
def createDescr():
    descr = InputFormDescr(title = 'Testing InputForm')
    descr.append({'name':'checkbutton',
                  'widgetType':Tkinter.Checkbutton,
                  'wcfg':{'text':'Checkbutton',
                          'variable':Tkinter.IntVar()},
                  'gridcfg':{'sticky':'w'}})
    

    descr.append({'name':'radioselect',
                'widgetType':Pmw.RadioSelect,
                'listtext':['rb1', 'rb2', 'rb3','rb4', 'rb5', 'rb6',
                            'rb7', 'rb8', 'rb9','rb10', 'rb11', 'rb12'],
                'wcfg':{'labelpos':'n',
                        'label_text':'Radiobuttons: ',
                        'orient':'vertical',
                        'buttontype':'radiobutton'},
                'gridcfg':{'sticky':'w'} })

    descr.append({'name':'radioselect2',
                'widgetType':Pmw.RadioSelect,
                'listtext':['rb1', 'rb2', 'rb3','rb4', 'rb5', 'rb6',
                            'rb7', 'rb8', 'rb9','rb10', 'rb11', 'rb12'],
                'wcfg':{'labelpos':'n',
                        'label_text':'Radiobuttons: ',
                        'orient':'horizontal',
                        'buttontype':'radiobutton'},
                'gridcfg':{'sticky':'w'} })
    entries = [('Chocolate',None), ('Vanilla', None), ('Strawberry', None),
               ('Coffee',None), ('Pistachio', None), ('Rasberry',None),
               ]

    descr.append({'name':'listchooser',
                  'widgetType':ListChooser,
                  'wcfg':{'entries':entries},
                  'gridcfg':{'sticky':'w'}
                  })
    
    return descr

def tearDown():
    # del master
    master.destroy()

def test_inputform_1():
    """
    Function to test the instanciation of an InputForm object with the
    following arguments:
    master = Tkinter.Tk()
    root   = None
    descr  = InputFormDescr containing a checkbutton and a radioselect
             widget... very simple.
    and the default value for the other argument of the constructor without
    specifying them
    """
    descr = createDescr()
    form = InputForm(master, root, descr)
    setWidget = {'checkbutton':1, 'radioselect':'rb9', 'radioselect2':'rb4',
                 'listchooser':'Pistachio'}

    value = form.testForm(setWidget=setWidget)

    assert value['checkbutton']==1, 'Expected 1, got %s'%(value['checkbutton'])
    assert value['radioselect']=='rb9'
    assert value['radioselect2']=='rb4'
    assert value['listchooser']==['Pistachio',]
    form.destroy()
    
def test_inputform_2():
    """
    
    Function to test the instanciation of an InputForm object with the
    following arguments:
    master = Tkinter.Tk()
    root   = None
    descr  = InputFormDescr containing a checkbutton and a radioselect
             widget... very simple.
    and the default value for the other argument of the constructor
    modal    = 1
    blocking = 0,
    defaultDirection = 'row',
    closeWithWindow = 1,
    onDestroy = None
    okCfg     = {'text':'OK'}
    cancelCfg = {'text':'Cancel'}
    initFunc  = None

    Focus on the scrolled frame options.
    """
    import Pmw
    descr = createDescr()
    form = InputForm(master, root, descr, modal=1, blocking=0,
                     defaultDirection='row', closeWithWindow=1,
                     onDestroy=None, okCfg={'text':'Validate'},
                     cancelCfg={'text':'Dismiss'},
                     scrolledFrame = 1,
                     width=400, height=500,
                     initFunc=None)
    
    assert isinstance(form.sf, Pmw.ScrolledFrame)
    setWidget = {'checkbutton':1, 'radioselect':'rb9', 'radioselect2':'rb4',
                 'listchooser':'Pistachio'}

    value = form.testForm(setWidget=setWidget)
    form.destroy() 


def test_inputform_3():
    """
    Function to test the instanciation of an InputForm object with the
    following argument, same options than before but this time the default
    arguments are specified. Change the text of the OK button to VALIDATE and
    the cancel button to DISMISS
    
    master = Tkinter.Tk()
    root   = None
    descr  = InputFormDescr containing a checkbutton and a radioselect
             widget... very simple.
    modal    = 1
    blocking = 0,
    defaultDirection = 'row',
    closeWithWindow = 1,
    onDestroy = None,
    sFrameCfg = {}
    okCfg     = {'text':'VALIDATE'}
    cancelCfg = {'text':'DISMISS'}
    initFunc  = None

    Focus on the ok and cancel button options
    """
    descr = createDescr()
    # I don't know how to make sure that the ok button is now labeled VALIDATE
    # and the cancel button 'DISMISS'.
    form = InputForm(master, root, descr, modal=1, blocking=0,
                     defaultDirection='row', closeWithWindow=1,
                     onDestroy=None, okCfg={'text':'VALIDATE'},
                     cancelCfg={'text':'DISMISS'},
                     initFunc=None)
    setWidget = {'checkbutton':1, 'radioselect':'rb9', 'radioselect2':'rb4',
                 'listchooser':'Pistachio'}

    value = form.testForm(setWidget=setWidget)
    #value = form.testForm()
    form.destroy() 
    
def test_inputform_okcfg():
    """
    Function to test the creation of an object InputForm with:
    master = Tkinter.Tk()
    root   = None
    descr  = InputFormDescr containing a checkbutton and a radioselect
             widget... very simple.
    okCfg     = {'text':'Validate', 'command':validate_cb}
    cancelCfg = {'text':'Cancel', 'command':dismiss_cb}
    The other argument keep their default values.
    
    Focus on the ok and cancel buttons option. Adding a function to be called
    after by ok_cb or cancel_cb. In this case the additional functions
    dismiss_cb and validate_cb don't take any arguments
    """
    global dismiss
    global validate
    dismiss = False
    validate = False
    def dismiss_cb():
        global dismiss
        dismiss = True

    def validate_cb():
        print "In validate_cb"
        global validate
        validate = True
        print validate
        
    descr = createDescr()
    form = InputForm(master, root, descr, modal=1, blocking=0,
                     okCfg={'text':'Validate',
                            'command':validate_cb},
                     cancelCfg={'text':'Dismiss',
                                'command':dismiss_cb})
    setWidget = {'checkbutton':1, 'radioselect':'rb9', 'radioselect2':'rb4',
                 'listchooser':'Pistachio'}
    value = form.testForm(setWidget=setWidget)
    print 'validate', validate
    print 'dismiss', dismiss
    assert validate
    assert not dismiss
    form.destroy() 

    
def test_inputform_initialize():
    """
    Function to test the creation of an object InputForm with:
    master = Tkinter.Tk()
    root   = None
    descr  = InputFormDescr containing a checkbutton and a radioselect
             widget... very simple.
    and the default value for the other argument of the constructor
    initFunc  = initialize()
    Focus on the initFunc
    """
    def initialize(form):
        w = form.descr.entryByName['listchooser']['widget']
        w.add(('Coconut', None))

    descr = createDescr()
    form = InputForm(master, root, descr, 
                     initFunc=initialize)
    setWidget = {'checkbutton':1, 'radioselect':'rb9', 'radioselect2':'rb4',
                 'listchooser':'Pistachio'}
    value = form.testForm(setWidget=setWidget)
    ent = form.descr.entryByName['listchooser']['widget'].entries
    assert 'Coconut' in map(lambda x: x[0], ent)

def test_inputform_groupwidgetsdefault():
    descr = InputFormDescr(title = 'Testing InputForm')
    descr.append({'name':'group',
                  'widgetType':Tkinter.Radiobutton,
                  'listtext':['rb1', 'rb2', 'rb3'],
                  'defaultValue':'rb3',
                  'gridcfg':{'sticky':'w'}})
    
    form = InputForm(master, root, descr)
    value = form.testForm()
    assert value['group']=='rb3'
    form.destroy() 

## def test_inputform_groupwidgets():
##     descr = InputFormDescr(title = 'Testing InputForm')
##     descr.append({'name':'group',
##                   'widgetType':Tkinter.Radiobutton,
##                   'listtext':['rb1', 'rb2', 'rb3'],
##                   'defaultValue':'rb3',
##                   'wcfg':{'variable':Tkinter.StringVar()},
##                   'gridcfg':{'sticky':'w'}})
    
##     form = InputForm(master, root, descr)
##     setWidget={'group':'rb1'}
##     value = form.testForm(setWidget=setWidget)
##     print value['group']
##     assert value['group']=='rb1'

def test_thumbwheel():
    def tw_cb(event=None):
        pass
    from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
    descr = InputFormDescr(title = 'Testing InputForm')
    descr.append({'name':'thumbwheel',
                  'widgetType':ThumbWheel,
                  'tooltip':"""Right click on the widget will display the control
                  panel and you can type a value manually""",
                  'defaultValue':100,
                  'wcfg':{'text':None,
                          'showLabel':1, 'width':100,
                          'min':0,
                          'lockBMin':1,
                          'lockBMax':1,
                          'lockBIncrement':1,
                          'value':40,
                          'oneTurn':1000,
                          'type':'int',
                          'increment':2,
                          'canvascfg':{'bg':'red'},
                          'wheelLabcfg1':{'font':(ensureFontCase('times'),14,'bold')},
                          'wheelLabcfg2':{'font':(ensureFontCase('times'),14,'bold')},
                          'callback':tw_cb,
                          'continuous':1, 'wheelPad':1, 'height':20},
                  'gridcfg':{'sticky':'e','row':-1}})
    
    form = InputForm(master, root, descr)
    value = form.testForm()
    assert value['thumbwheel']==100
    form.destroy() 

def test_scrolledText():
    descr = InputFormDescr(title = "Testing ScrolledText")
    descr.append({'widgetType':Pmw.ScrolledText,
                  'name':'sText',
                  'defaultValue':"""DEFAULT TEXT""",
                  'wcfg':{'labelpos':'n',
                          'label_text':'ScrolledText with headers',
                          'usehullsize': 1,
                          'hull_width': 400,
                          'hull_height': 300,
                          'text_wrap':'none',
                          'text_padx': 4,
                          'text_pady': 4,
                          },
                  'gridcfg':{'sticky':'wens'}})
    
    form = InputForm(master, root, descr, modal=0, blocking=0)
    values = form.testForm()
    if not values['sText'] == 'DEFAULT TEXT\n':
        raise RuntimeError
    form.destroy() 

                
def test_inputForm_notebook():
    descr = InputFormDescr(title = 'Testing InputForm')
    descr.append({'widgetType':Pmw.NoteBook,
                  'name':'notebook',
                  'container':{'Page1':"w.page('Page1')",
                               'Page2':"w.page('Page2')",
                               'Page3':"w.page('Page3')"},
                  'wcfg':{'borderwidth':3},
                  'componentcfg':[{'name':'Page1', 'cfg':{}},
                                  {'name':'Page2', 'cfg':{}},
                                  {'name':'Page3', 'cfg':{}}],
                  'gridcfg':{'sticky':'wnse'}
                  })

    entries = [('Chocolate',None), ('Vanilla', None), ('Strawberry', None),
               ('Coffee',None), ('Pistachio', None), ('Rasberry',None),
               ]

    descr.append({'name':'listchooser',
                  'parent':'Page1',
                  'widgetType':ListChooser,
                  'defaultValue':'Chocolate',
                  'wcfg':{'entries':entries},
                  'gridcfg':{'sticky':'w'}
                  })
    
    descr.append({'name':'radioselect2',
                  'widgetType':Pmw.RadioSelect,
                  'parent':'Page2',
                  'listtext':['rb1', 'rb2', 'rb3','rb4', 'rb5', 'rb6',
                              'rb7', 'rb8', 'rb9','rb10', 'rb11', 'rb12'],
                  'wcfg':{'labelpos':'n',
                          'label_text':'Radiobuttons: ',
                          'orient':'horizontal',
                          'buttontype':'radiobutton'},
                  'gridcfg':{'sticky':'w'} })

    descr.append({'name':'radioselect',
                  'widgetType':Pmw.RadioSelect,
                  'parent':'Page3',
                  'defaultValue':'rb5',
                  'listtext':['rb1', 'rb2', 'rb3','rb4', 'rb5', 'rb6',
                              'rb7', 'rb8', 'rb9','rb10', 'rb11', 'rb12'],
                  'wcfg':{'labelpos':'n',
                          'label_text':'Radiobuttons: ',
                          'orient':'vertical',
                          'buttontype':'radiobutton'},
                  'gridcfg':{'sticky':'w'} })
    form = InputForm(master, root, descr, modal=0, blocking=0)
    values = form.testForm(container='Page3')
    assert values['radioselect']=='rb5'
    form.destroy() 
    
harness = testplus.TestHarness( __name__,
                                connect = setUp,
                                funs = testplus.testcollect( globals()),
                                #funs = [test_inputform_4,] ,
                                disconnect = tearDown
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))
